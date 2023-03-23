import pickle
from pathogenprofiler import cmd_out,debug,infolog,errlog
import os
import json
from .output import write_outputs
from copy import copy
import filelock
from time import time
import sqlite3
from tqdm import tqdm

def extract_variant_set(vcf_file, exclude_bed, min_cov=10, min_freq=0.8):
    ref_diffs = set()
    missing = set()
    for l in cmd_out(f"bcftools view -V indels -T ^{exclude_bed} {vcf_file} | bcftools query -f '%POS[\t%GT:%AD]\n'"):
        if l[0]=="#": continue
        row = l.strip().split()
        pos = int(row[0])
        gt,ad = row[1].split(":")
        if ad==".": # delly
            continue
        if gt==".": 
            missing.add(pos)
            continue
        ad = [int(x) for x in ad.split(",")]
        if sum(ad)<=min_cov:
            missing.add(pos)
            continue
        adf = sorted([float(x/sum(ad)) for x in ad])
        if adf[-1]<min_freq:
            missing.add(pos)
            continue
        if gt=="1/1":
            ref_diffs.add(int(pos))

    return ref_diffs,missing

class DB:
    def __init__(self, filename):
        self.filename = filename
        self.conn = sqlite3.connect(filename)
        self.c = self.conn.cursor()
        self.c.execute('''CREATE TABLE IF NOT EXISTS variants (sample text, lineage text, diffs binary, missing binary)''')
    def store(self,json_results, vcf_file, exclude_bed, min_cov=10, min_freq=0.8):
        sample_name = json_results["id"]
        diffs,missing = extract_variant_set(vcf_file, exclude_bed, min_cov, min_freq)
        res = self.c.execute("SELECT sample FROM variants WHERE sample=?",(sample_name,)).fetchone()
        if res:
            self.c.execute("UPDATE variants SET lineage=?, diffs=?, missing=? WHERE sample=?",(json_results['sublin'], pickle.dumps(diffs), pickle.dumps(missing), sample_name))
        else:
            self.c.execute("INSERT INTO variants VALUES (?,?,?,?)",(sample_name, json_results['sublin'], pickle.dumps(diffs), pickle.dumps(missing)))
        self.conn.commit()
        self.diffs = diffs
        self.missing = missing
    def search(self,json_results,vcf_file, exclude_bed, cutoff = 20, min_cov=10, min_freq=0.8):
        debug("Searching for close samples in %s" % self.filename)
        self.c.execute("SELECT sample, diffs, missing FROM variants WHERE lineage=?",(json_results['sublin'],))
        self.diffs,self.missing = extract_variant_set(vcf_file, exclude_bed, min_cov, min_freq)
        sample_dists = []
        for s,d,m in tqdm(self.c.fetchall()):
            dist = self.diffs.symmetric_difference(pickle.loads(d))
            dist -= self.missing
            dist -= pickle.loads(m) 
            if (ld:=len(dist))<cutoff:
                snp_union = self.diffs.union(pickle.loads(d))
                missing_union = self.missing.union(pickle.loads(m))
                percent_missing = len(snp_union.intersection(missing_union))/len(snp_union)
                sample_dists.append({
                    "sample":s,
                    "distance":ld,
                    "diffs":list(dist),
                    "percent_missing": percent_missing
                })
        return sample_dists

def read_json(filename):
    debug("Reading %s" % filename)
    lock = filelock.FileLock(filename + ".lock")
    with lock:
        data = json.load(open(filename))
    debug("Finished reading %s" % filename)
    return data

def sample_in_json(sample,result_file):
    results = read_json(result_file)
    if sample in [d['sample'] for d in results['close_samples']]:
        return True
    else:
        return False

def run_snp_dists(args,results):
    infolog("\nCalculating SNP distances to look for closely related samples")
    infolog("-------------------------------------------------------------")
    t1 = time()
    if args.vcf:
        wg_vcf = args.vcf
    else:
        wg_vcf = args.files_prefix + ".vcf.gz"
    if args.snp_diff_db:
        dbname = args.snp_diff_db
    else:
        dbname = f'{args.dir}/results/snp_diffs.db'
    db = DB(dbname)
    if not args.snp_diff_no_store:
        db.store(results,wg_vcf,args.conf['bedmask'])
    results["close_samples"] = db.search(results,wg_vcf,args.conf['bedmask'],args.snp_dist)
    results["close_samples"] = [d for d in results["close_samples"] if d["sample"]!=results["id"]]
    t2 = time()
    dt = int(t2-t1)
    infolog(f"\nFound {len(results['close_samples'])} close samples in {dt} seconds")
    infolog("-------------------------------------------------------------\n")
    # if args.nj:
        # results['_tree'] = make_nj_tree(args,results)

def update_neighbour_snp_dist_output(args,results):
    for s in results['close_samples']:
        infolog("Close sample found: %s (%s). Updating result files" % (s['sample'],s['distance']))
        f = os.path.join(args.dir,"results",f"{s['sample']}.results.json")
        if not os.path.exists(f):
            continue
        if not sample_in_json(args.prefix,f):
            lock = filelock.FileLock(f + ".lock")
            with lock:
                debug("Acquiring lock for %s" % f)
                data = json.load(open(f))
                data['close_samples'].append({
                    "sample":args.prefix,
                    "distance":s['distance']
                })
                temp_args = copy(args)
                temp_args.prefix = s['sample']
                # if args.nj:
                    # results['_tree'] = make_nj_tree(temp_args,data)
                write_outputs(temp_args,data,template_file=args.text_template)
                debug("Finished with lock for %s" % f)

# def make_nj_tree(args,results):
    
#     neighbours = [s['sample'] for s in results['close_samples']]
#     if len(neighbours)<2:
#         return None
#     all_samps = neighbours + [results['id']]
#     if len(neighbours)==0: return
#     dist_dict = {}
#     for d in results['close_samples']:
#         dist_dict[(results['id'],d['sample'])] = d['distance']
#         dist_dict[(d['sample'],results['id'])] = d['distance']
#     for si in neighbours:
#         data = read_json(os.path.join(args.dir,"results",f"{si}.results.json"))
#         for d in data['close_samples']:
#             sj = d['sample']
#             if sj not in all_samps: continue
#             dist_dict[(si,sj)] = d['distance']
#             dist_dict[(sj,si)] = d['distance']
#     dists = []

#     for si in all_samps:
#         row = []
#         for sj in all_samps:
#             if si==sj: 
#                 row.append(0)
#             elif (si,sj) in dist_dict:
#                 row.append(dist_dict[(si,sj)])
#             else:
#                 si_set_file = os.path.join(args.dir,"results",f"{si}.non_ref.pkl")
#                 sj_set_file = os.path.join(args.dir,"results",f"{sj}.non_ref.pkl")
#                 si_set = variant_set(si_set_file)
#                 d = si_set.get_snp_dist(sj_set_file)
#                 row.append(d)
#         dists.append(row)
#     from skbio import DistanceMatrix
#     from skbio.tree import nj
#     dm = DistanceMatrix(dists, all_samps)
#     tree = nj(dm)
#     tree = tree.root_at_midpoint()
#     return tree


