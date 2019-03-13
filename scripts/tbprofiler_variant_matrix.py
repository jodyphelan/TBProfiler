#! /usr/bin/env python
from __future__ import division
import pathogenprofiler as pp
import sys
import json
import csv
import argparse
from collections import defaultdict
import math
from tqdm import tqdm

def main(args):
	sample_file = args.samples
	dst_file = args.dst

	dst = json.load(open(dst_file))

	samples = [x.rstrip() for x in open(sample_file).readlines()]
	ext = ".results.json"
	var = set()
	variants = defaultdict(lambda:defaultdict(dict))
	for s in tqdm(samples):
		res_file = "%s%s" % (s,ext)
		if pp.nofile(res_file):
			pp.log("Warning: %s does not exist!" % res_file)
			continue
		res = json.load(open(res_file))
		for v in res["dr_variants"]:
			var.add((v["locus_tag"],v["change"]))
			variants[v["locus_tag"]][v["change"]][s] = "1"
		for v in res["other_variants"]:
			var.add((v["locus_tag"],v["change"]))
			variants[v["locus_tag"]][v["change"]][s] = "1"
	outfile = "%s.matrix.txt" % (args.prefix)
	O = open(outfile,"w")
	O.write("gene\tvariant\t%s\n" % ("\t".join(samples)))
	for gene,change in tqdm(var):
		for s in samples:
			if s not in variants[gene][change]:
				variants[gene][change][s] = "0"
		O.write("%s\t%s\t%s\n" % (gene,change,"\t".join([variants[gene][change][s] for s in samples])))
	O.close()

	drug_loci = pp.load_bed(args.bed,[6],4) #{'rrl': ['linezolid']}
	for d in drug_loci:
		drug_loci[d] = [x.lower() for x in drug_loci[d][0].split(";")]

	for gene,change in tqdm(var):
		for drug in drug_loci[gene]:
			if drug not in dst[s]: continue
			a,b,c,d = 0.5,0.5,0.5,0.5
			for s in samples:
				if dst[s][drug]=="1" and variants[gene][change][s]=="1": a+=1 # 1 R
				if dst[s][drug]=="0" and variants[gene][change][s]=="1": b+=1 # 1 S
				if dst[s][drug]=="1" and variants[gene][change][s]=="0": c+=1 # 0 R
				if dst[s][drug]=="0" and variants[gene][change][s]=="0": d+=1 # 0 S
			OR = (a/b)/(c/d)
			SE = math.sqrt((1/a)+(1/b)+(1/c)+(1/d))
			Q = math.exp(1.96*SE)
			L_CI = OR/Q
			H_CI = OR*Q
			print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (gene,change,drug,a,b,c,d,OR,L_CI,H_CI))
#	json.dump(variants,open("debug.json","w"))



parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_sub = subparsers.add_parser('matrix', help='Run whole pipeline', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('samples',help='NGS Platform')
parser_sub.add_argument('dst',help='NGS Platform')
parser_sub.add_argument('prefix',help='NGS Platform')
parser_sub.add_argument('bed',help='NGS Platform')
parser_sub.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
