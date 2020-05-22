#! /usr/bin/env python
import json
from collections import defaultdict
import re
import argparse
import os
from tqdm import tqdm
import sys

def main(args):
	change_field = "change" if args.variant_format=="hgvs" else "_internal_change"
	conf = json.load(open(sys.prefix+"/share/tbprofiler/%s.config.json" % args.db))
	drug2genes = defaultdict(set)
	for l in open(conf["bed"]):
		row = l.rstrip().split()
		for d in row[5].split(","):
			drug2genes[d].add(row[3])
	if args.samples:
		samples = [x.rstrip() for x in open(args.samples).readlines()]
	else:
		samples = [x.replace(".results.json","") for x in os.listdir("results/") if x[-13:]==".results.json"]
	variants = defaultdict(lambda:defaultdict(list))
	meta = json.load(open(args.meta))
	mutation_types = defaultdict(dict)
	for s in tqdm(samples):
		tmp = json.load(open("%s/%s.results.json" % (args.dir,s)))
		for var in tmp["dr_variants"]:
			variants[var["locus_tag"]][var[change_field]].append(s)
			mutation_types[var["locus_tag"]][var[change_field]] = var["type"]
		for var in tmp["other_variants"]:
			variants[var["locus_tag"]][var[change_field]].append(s)
			mutation_types[var["locus_tag"]][var[change_field]] = var["type"]
	print("drug\tgene\tmutation\ttype\tmutation_whole_dataset\tnum_samples_whole_dataset\tmutation_dst\tno_mutation_dst\tOR")
	for drug in sorted(drug2genes):
		if drug not in meta[samples[0]]: continue
		for gene in sorted(drug2genes[drug]):
			if gene not in variants: continue
			for mutation in tqdm(variants[gene]):
				t = [
					[0.5,0.5],
					[0.5,0.5]
					  ]
				for s in samples:
					if s not in meta: continue
					if meta[s][drug]=="1" and s in variants[gene][mutation]: 		t[0][0]+=1
					if meta[s][drug]=="0" and s in variants[gene][mutation]: 		t[0][1]+=1
					if meta[s][drug]=="1" and s not in variants[gene][mutation]: 	t[1][0]+=1
					if meta[s][drug]=="0" and s not in variants[gene][mutation]: 	t[1][1]+=1
				OR = (t[0][0]/t[0][1])/(t[1][0]/t[1][1])
				if sum(t[0])-1==0:
					OR = "NA"
				print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (drug,gene,mutation,mutation_types[gene][mutation],len(variants[gene][mutation]),len(samples),int(sum(t[0])-1),int(sum(t[1])-1),OR))
parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('meta',type=str,help="Meta data")
parser.add_argument('--samples',type=str,help='NGS Platform')
parser.add_argument('--dir',default="results/",type=str,help='NGS Platform')
parser.add_argument('--db',default="tbdb",type=str,help='NGS Platform')
parser.add_argument('--variant-format',default="hgvs",choices=["hgvs","bcftools"],type=str,help='NGS Platform')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
