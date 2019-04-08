#! /usr/bin/env python
import json
import sys
from collections import defaultdict
import re
import argparse

def main(args):
	drugs = [
	'isoniazid', 'rifampicin', 'ethambutol', 'pyrazinamide', 'streptomycin',
	'fluoroquinolones', 'amikacin', 'capreomycin', 'kanamycin',
	'cycloserine',  'ethionamide', 'clofazimine', 'para-aminosalicylic_acid',
	'delamanid', 'bedaquiline', 'linezolid', ]
	if args.meta:
		meta = json.load(open(args.meta))
		col,val = args.meta_col.split(",")
		meta_samples = [s for s in meta if mets[s][col]==val]
	variants = defaultdict(lambda:defaultdict(int))
	data = json.load(open(sys.argv[1]))
	for s in data:
		if args.meta and s not in meta_samples: continue
		for d in drugs:
			if data[s][d]=="-": continue
			muts = [x.strip() for x in data[s][d].split(",")]
			for m in muts:
				variants[d][m]+=1

	for d in drugs:
		for m in variants[d]:
			re_obj = re.search("([a-zA-Z0-9]+)_(.*)",m)
			gene = re_obj.group(1)
			var = re_obj.group(2)
			print("%s\t%s\t%s\t%s" % (d.capitalize(),gene,var,variants[d][m]))

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('json',help='NGS Platform')
parser.add_argument('--meta',type=str)
parser.add_argument('--meta-col',type=str)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
