#! /usr/bin/env python
import os
import sys
import argparse
import json
from collections import defaultdict
from tqdm import tqdm

def main(args):
	if args.samples:
		samples = [x.rstrip() for x in open(args.sample_file).readlines()]
	else:
		samples = [x.replace(".results.json","") for x in os.listdir("%s/"%args.dir) if x[-13:]==".results.json"]
	mutations = defaultdict(set)
	for s in tqdm(samples):
		tmp = json.load(open("%s/%s.results.json" % (args.dir,s)))
		for var in tmp["dr_variants"]:
			mutations[(var["gene"],var["change"])].add(s)
		for var in tmp["other_variants"]:
			mutations[(var["gene"],var["change"])].add(s)
	O = open(args.out,"w") if args.out else sys.stdout
	O.write("gene\tchange\t%s\n" % ("\t".join(samples)))
	for gene,change in tqdm(mutations):
		O.write("%s\t%s\t%s\n" % (gene,change,"\t".join(["1" if s in mutations[(gene,change)] else "0" for s in samples])))
	O.close()


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--samples',type=str,help='NGS Platform')
parser.add_argument('--dir',default="results/",type=str,help='NGS Platform')
parser.add_argument('--out',default="results/",type=str,help='NGS Platform')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
