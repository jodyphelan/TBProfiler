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
	print("sample,%s" % (",".join(["mutations_%s" % x for x in sorted(drug2genes)])))
	for s in tqdm(samples):
		tmp = json.load(open("%s/%s.results.json" % (args.dir,s)))
		sample_mutations = defaultdict(list)
		tmp_store = set([json.dumps(x) for x in tmp["dr_variants"]])
		tmp["dr_variants"] = []
		tmp["dr_variants"] = [json.loads(x) for x in tmp_store]
		for var in tmp["dr_variants"]:
			sample_mutations[var["drug"]].append(var["gene"]+"_"+var[change_field])
		print("%s,%s" % (s,",".join(["; ".join(sample_mutations[d]) for d in sorted(drug2genes)])))

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--samples',type=str,help='NGS Platform')
parser.add_argument('--dir',default="results/",type=str,help='NGS Platform')
parser.add_argument('--db',default="tbdb",type=str,help='NGS Platform')
parser.add_argument('--variant-format',default="hgvs",choices=["hgvs","bcftools"],type=str,help='NGS Platform')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
