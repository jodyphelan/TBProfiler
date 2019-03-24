import sys
import argparse
import os
import json
from collections import defaultdict
from tqdm import tqdm
def main(args):
	freqs = defaultdict(lambda:defaultdict(list))

	samples = [x.replace(".results.json","") for x in os.listdir(args.dir) if x[-13:]==".results.json"]
	for s in tqdm(samples):
		tmp = set()
		results = json.load(open("%s/%s.results.json" % (args.dir,s)))
		for var in results["dr_variants"]:
			tmp.add((var["gene"],var["change"],var["freq"]))
			if var["gene"]=="Rv0678":
				print("%s\t%s\t%s\n" % (s,var["change"],var["freq"]))
		for gene,change,freq in tmp:
			freqs[gene][change].append(freq)

	print("Gene\tHomozygous\tHeterozygous")
	for gene in freqs:
		hom = 0
		het = 0
		for mut in freqs[gene]:
			for f in freqs[gene][mut]:
				if f<args.cutoff:
					het+=1
				else:
					hom+=1
		print("%s\t%s\t%s" % (gene,hom,het))

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--dir','-d',default="results/",type=str,help='NGS Platform')
parser.add_argument('--cutoff','-c',default=0.95,type=float,help='NGS Platform')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
