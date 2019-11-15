#! /usr/bin/env python
import os
import sys
import argparse
import json
from collections import defaultdict
from tqdm import tqdm
import pathogenprofiler as pp
def main(args):
	total_samples = set()
	sample_sets = {}
	for f in args.samples.split(","):
		sample_sets[f] = [s.rstrip() for s in open(f).readlines()]
		for s in sample_sets[f]:
			total_samples.add(s)

	mutations = defaultdict(set)
	columns = ["gene","change"] + (args.columns.split(",") if args.columns else [])

	for s in tqdm(total_samples):
		if pp.nofile("%s/%s.results.json" % (args.dir,s)): continue
		tmp = json.load(open("%s/%s.results.json" % (args.dir,s)))
		for var in tmp["dr_variants"]:
			tmp_var = {x:var[x] for x in columns}
			mutations[json.dumps(tmp_var)].add(s)
		if args.non_dr:
			for var in tmp["other_variants"]:
				tmp_var = var
				tmp_var = {x:var[x] for x in columns}
				mutations[json.dumps(tmp_var)].add(s)


	if args.summary:
		O = open(args.summary,"w")
		O.write("%s\n" % ("\t".join(["Gene","Mutation"]+columns+["%s_num\t%s_pct" % (x,x) if args.pct else x for x in list(sample_sets)])))
		for var_string in mutations:
			if "gene_name" not in var: var["gene_name"] = var["gene"] ######Fix for large deletions not haveing this key

			tmp_freqs = []
			for f in sample_sets:
				num = len([s for s in sample_sets[f] if s in mutations[var_string]])
				if args.pct:
					pct = num/len(sample_sets[f])*100
					tmp_freqs.append("%s\t%.2f" % (num,pct))
				else:
					tmp_freqs.append("%s" % num)
			O.write("%s\t%s\t%s\t%s\n" % (var["gene_name"],var["change"],"\t".join([str(var[x]) for x in columns]),"\t".join(tmp_freqs)))
		O.close()

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('samples',type=str,help='NGS Platform')
parser.add_argument('--dir',default="results",type=str,help='NGS Platform')
parser.add_argument('--summary',type=str,help='NGS Platform')
parser.add_argument('--pct',action="store_true",help='NGS Platform')
parser.add_argument('--columns',type=str,help='NGS Platform')
parser.add_argument('--non-dr',action="store_true",help='NGS Platform')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
