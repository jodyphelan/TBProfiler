#! /usr/bin/env python
import json
import sys
from collections import defaultdict
import re
import argparse
from collections import defaultdict

def main(args):
	drugs = [
	'isoniazid', 'rifampicin', 'ethambutol', 'pyrazinamide', 'streptomycin',
	'fluoroquinolones', 'amikacin', 'capreomycin', 'kanamycin',
	'cycloserine',  'ethionamide', 'clofazimine', 'para-aminosalicylic_acid',
	'delamanid', 'bedaquiline', 'linezolid', ]
	if args.meta:
		meta = json.load(open(args.meta))
		samples2meta = {s:meta[s][args.meta_col] for s in meta}
		meta_cats = {x:0 for x in samples2meta}
		meta_cats["NA"] = 0
	variants = defaultdict(lambda:defaultdict(int))
	data = json.load(open(sys.argv[1]))
	for s in data:
#		if args.meta and s not in meta_samples: continue
		for d in drugs:
			if data[s][d]=="-": continue
			muts = [x.strip() for x in data[s][d].split(",")]
			for m in muts:
				variants[d][m]+=1
				if args.meta:
					if s in samples2meta:
						meta_cats[samples2meta[s]]+=1
					else:
						meta_cats["NA"]+=1

	for d in drugs:
		for m in variants[d]:
			re_obj = re.search("([a-zA-Z0-9]+)_(.*)",m)
			gene = re_obj.group(1)
			var = re_obj.group(2)
			sys.stdout.write("%s\t%s\t%s\t%s" % (d.capitalize(),gene,var,variants[d][m]))
			if args.meta:
				for cat in meta_cats:
					num = meta_cats[cat]
					tot_num = len([x for x in samples2meta if samples2meta[x]==cat])
					pct = num/tot_num*100
					sys.stdout.write("\t%s/%s (%.2f)" % (num,tot_num,pct))
			sys.stdout.write("\n")

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('json',help='NGS Platform')
parser.add_argument('--meta',type=str)
parser.add_argument('--meta-col',type=str)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
