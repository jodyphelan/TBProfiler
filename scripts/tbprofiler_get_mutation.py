#! /usr/bin/env python
import json
import sys
import argparse
from tqdm import tqdm
import os
def main(args):#locus_tag,mutation,resultdir=".",sample_file=None):
	if args.sample_file:
		samples = [x.rstrip() for x in open(args.sample_file).readlines()]
	else:
		samples = [x.replace(".results.json","") for x in os.listdir("%s/"%args.dir) if x[-13:]==".results.json"]
	positives = []
	for s in tqdm(samples):
		tmp = json.load(open("%s/%s.results.json" % (args.dir,s)))
		for var in tmp["dr_variants"]:
			if var["gene"]==args.gene and var["change"]==args.mutation:
				positives.append(s)
		for var in tmp["other_variants"]:
			if var["gene"]==args.gene and var["change"]==args.mutation:
				positives.append(s)

	if args.txt:
		O = open("%s_%s.samples.txt" % (args.gene,args.mutation.replace(">","_")),"w")
		O.write("%s\n" %("\n".join(positives)))
		O.close()
	if args.itol:
		O = open("%s_%s.itol.txt" % (args.gene,args.mutation.replace(">","_")),"w")
		O.write("""DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL	samples
COLOR	#ff0000

LEGEND_TITLE	samples
LEGEND_SHAPES	1
LEGEND_COLORS	black
LEGEND_LABELS	Sample

DATA
""")
		for x in positives:
			O.write("%s\tblack\n" % x.rstrip())
		O.close()
	if args.lineage:
		lineages = {}
		for l in open(args.lineage):
			row = l.rstrip().split()
			print(row)
			lineages[row[0]] = row[1]
		for s in positives:
			sys.stdout.write("%s\t%s\n" % (s,lineages[s]))
	if not args.lineage and not args.txt and not args.itol:
		sys.stdout.write("%s\n" % "\n".join(positives))

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('gene',help='NGS Platform')
parser.add_argument('mutation',help='NGS Platform')
parser.add_argument('--sample_file','-s',default=None,type=str,help='NGS Platform')
parser.add_argument('--dir','-d',default="results/",type=str,help='NGS Platform')
parser.add_argument('--lineage','-l',default=None,type=str,help='NGS Platform')
parser.add_argument('--txt',action='store_true')
parser.add_argument('--itol',action='store_true')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
