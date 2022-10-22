import quickle
from pathogenprofiler import cmd_out,debug,infolog
import os
import json
from .output import write_outputs
from copy import copy
import filelock
from time import time

def write_variant_set(vcf_file, prefix, exclude_bed, min_cov=10, min_freq=0.8):
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

    filename = prefix + ".non_ref.qkl"
    open(filename,"wb").write(quickle.dumps((ref_diffs,missing)))
    return variant_set(filename)

class variant_set:
    def __init__(self, filename):
        self.filename = filename
        self.diffs,self.missing = quickle.loads(open(filename,"rb").read())
    
    def get_snp_dist(self, set_file):
        other_diffs,other_missing = quickle.loads(open(set_file,"rb").read())
        pairwise_dists = self.diffs.symmetric_difference(other_diffs)
        pairwise_dists -= self.missing
        pairwise_dists -= other_missing
        return len(pairwise_dists)
    def get_close_samples(self,dir,cutoff=20):
        self.sample_dists = []
        directory_files = [f for f in os.listdir(dir) if f.endswith(".qkl")]
        infolog(f"Searching across {len(directory_files)} files")
        for f in directory_files:
            dist = self.get_snp_dist(os.path.join(dir,f))
            self.sample_dists.append({
                "sample":f.replace(".non_ref.qkl",""),
                "distance":dist
            })
        return [x for x in self.sample_dists if x["distance"]<=cutoff]

def sample_in_json(sample,result_file):
    results = json.load(open(result_file))
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
    var_set = write_variant_set(wg_vcf,args.files_prefix,exclude_bed=args.conf['bedmask'])
    results["close_samples"] = var_set.get_close_samples(os.path.join(args.dir,"results"),cutoff=args.snp_dist)
    results["close_samples"] = [d for d in results["close_samples"] if d["sample"]!=results["id"]]
    i=0
    for s in results['close_samples']:
        infolog("Close sample found: %s (%s). Updating result files" % (s['sample'],s['distance']))
        f = os.path.join(args.dir,"results",f"{s['sample']}.results.json")
        if not sample_in_json(args.prefix,f):
            lock = filelock.FileLock(f + ".lock")
            with lock:
                data = json.load(open(f))
                data['close_samples'].append({
                    "sample":args.prefix,
                    "distance":s['distance']
                })
                temp_args = copy(args)
                temp_args.prefix = s['sample']
                write_outputs(temp_args,data,template_file=args.output_template)
    t2 = time()
    dt = int(t2-t1)
    infolog(f"\nFound {len(results['close_samples'])} close samples in {dt} seconds")
    infolog("-------------------------------------------------------------\n")


def make_nj_tree(args,results):
    from skbio import DistanceMatrix
    from skbio.tree import nj
    neighbours = [s['sample'] for s in results['close_samples']]
    all_samps = neighbours + [results['id']]
    if len(neighbours)==0: return
    dist_dict = {}
    for d in results['close_samples']:
        dist_dict[(results['id'],d['sample'])] = d['distance']
        dist_dict[(d['sample'],results['id'])] = d['distance']
    for si in neighbours:
        data = json.load(open(os.path.join(args.dir,"results",f"{si}.results.json")))
        for d in data['close_samples']:
            sj = d['sample']
            if sj not in all_samps: continue
            dist_dict[(si,sj)] = d['distance']
            dist_dict[(sj,si)] = d['distance']
    dists = []

    for si in all_samps:
        row = []
        for sj in all_samps:
            if si==sj: 
                row.append(0)
            elif (si,sj) in dist_dict:
                row.append(dist_dict[(si,sj)])
            else:
                si_set_file = os.path.join(args.dir,"results",f"{si}.non_ref.qkl")
                sj_set_file = os.path.join(args.dir,"results",f"{si}.non_ref.qkl")
                si_set = variant_set(si_set_file)
                d = si_set.get_snp_dist(sj_set_file)
                row.append(d)
        dists.append(row)

    dm = DistanceMatrix(dists, all_samps)
    tree = nj(dm)
    tree = tree.root_at_midpoint()
    return tree


