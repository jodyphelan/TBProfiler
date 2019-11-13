#! /usr/bin/env python
import json
from collections import defaultdict
import argparse
import os
from tqdm import tqdm
import sys
import pathogenprofiler as pp
import csv

try:
    sys.base_prefix
except:
    sys.base_prefix = getattr(sys, 'base_prefix', getattr(sys, 'real_prefix', sys.prefix))


def get_conf_dict(library_prefix):
    files = {"gff":".gff","ref":".fasta","ann":".ann.txt","barcode":".barcode.bed","bed":".bed","json_db":".dr.json","version":".version.json"}
    conf = {}
    for key in files:
        sys.stderr.write("Using %s file: %s\n" % (key,library_prefix+files[key]))
        conf[key] = pp.filecheck(library_prefix+files[key])
    return conf



def main(args):
    conf = conf = get_conf_dict(sys.base_prefix + "/share/tbprofiler/%s" % args.db)
    json_db = json.load(open(conf["json_db"]))
    drug2genes = defaultdict(set)
    gene2drugs = defaultdict(set)
    gene2lt = {}
    for l in open(conf["bed"]):
        row = l.rstrip().split()
        for d in row[5].split(","):
            drug2genes[d].add(row[3])
            gene2drugs[row[3]].add(d)
            gene2lt[row[3]] = row[3]
            gene2lt[row[4]] = row[3]

    mutations = set()
    for l in open(args.mutations):
        row = l.strip().split()
        mutations.add((gene2lt[row[0]], row[1]))

    meta = {}
    reader = csv.DictReader(open(args.meta))
    drug_resistant_isolates = {d:set() for d in drug2genes if d in reader.fieldnames}
    for row in reader:
        meta[row["id"]] = row
        for drug in drug_resistant_isolates:
            if row[drug]=="1":
                drug_resistant_isolates[drug].add(row["id"])
    pp.log(f"Analysing {len(drug_resistant_isolates)} drugs")
    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines() if x.rstrip() in meta]
    else:
        samples = [x.replace(".results.json","") for x in os.listdir("results/") if x[-13:]==".results.json" if x.replace(".results.json","") in meta]

    variants = {x:set() for x in mutations}
    hgvs2bcftools = {}
    for s in tqdm(samples):
        tmp = json.load(open(f"{args.dir}/{s}.results.json"))
        for var in tmp["dr_variants"]+tmp["other_variants"]:
            if (var["locus_tag"],var["change"]) in mutations:
                hgvs2bcftools[var["change"]] = var["_internal_change"]
                variants[(var["locus_tag"],var["change"])].add(s)
    total_sample_n = len(samples)
    pp.log(f"Found {total_sample_n} samples in meta list with result files")
    pp.log("-"*40)
    print("Gene,Mutation,Drug resistance association,Total frequency (percentage),Associated drugs,Drug resitant frequency (percentage)")
    for gene,mut in variants:
        if (gene,mut) not in variants:
            continue
        total_freq = len(variants[(gene,mut)])
        total_pct = total_freq/total_sample_n*100
        dr_associated = "Not associated"
        if gene in json_db and hgvs2bcftools[mut] in json_db[gene]:
            drugs = list(json_db[gene][hgvs2bcftools[mut]]["drugs"].keys())
            dr_associated = "Associated"
        else:
            drugs = gene2drugs[gene]
        dr_freqs = []
        dr_pcts = []
        for drug in drugs:
            dr_freq = len(variants[(gene,mut)].intersection(drug_resistant_isolates[drug]))
            dr_pct = dr_freq/len(drug_resistant_isolates[drug])*100
            dr_freqs.append(dr_freq)
            dr_pcts.append(dr_pct)

        zipped_list = ["%s (%.2f%%)" % (x,y) for x,y in zip(dr_freqs,dr_pcts)]
        print("%s,%s,%s,%s (%.2f),%s,%s" % (gene,mut,dr_associated,total_freq,total_pct,';'.join(drugs),';'.join(zipped_list)))
        pp.log("-"*40)

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--mutations',type=str,help='NGS Platform',required=True)
parser.add_argument('--meta',default="tbdb",type=str,help='NGS Platform',required=True)
parser.add_argument('--samples',type=str,help='NGS Platform')
parser.add_argument('--dir',default="results/",type=str,help='NGS Platform')
parser.add_argument('--db',default="tbdb",type=str,help='NGS Platform')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
