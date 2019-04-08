import sys
import argparse
import csv
from collections import defaultdict
import re
import json
from tqdm import tqdm
def main(args):

	drugs = [
		'isoniazid', 'rifampicin', 'ethambutol', 'pyrazinamide', 'streptomycin',
		'fluoroquinolones', 'amikacin', 'capreomycin', 'kanamycin',
		 'cycloserine',  'ethionamide', 'clofazimine', 'para-aminosalicylic_acid',
		 'delamanid', 'bedaquiline', 'linezolid', ]

	variants = defaultdict(lambda:defaultdict(int))
	data = json.load(open(args.json))
	for s in tqdm(data):
		for d in drugs:
			if data[s][d]=="-": continue
			muts = [x.strip() for x in data[s][d].split(",")]
			for m in muts:
				re_obj = re.search("([a-zA-Z0-9]+)_(.*)",m)
				gene = re_obj.group(1)
				var = re_obj.group(2)
				variants[gene][var]+=1
	IN = csv.DictReader(open(args.csv))
	OUT = open(args.out,"w")
	outwriter = csv.DictWriter(OUT,fieldnames = IN.fieldnames+["NGS"])
	for row in IN:
		tmp = row
		if row["Gene"] in variants and row["Mutation"] in variants[row["Gene"]]:
			tmp["NGS"] = "Yes" if args.binary else variants[row["Gene"]][row["Mutation"]]
		else:
			tmp["NGS"] = "No" if args.binary else 0
		outwriter.writerow(tmp)
	OUT.close()


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('csv',help='NGS Platform')
parser.add_argument('json',help='NGS Platform')
parser.add_argument('out',help='NGS Platform')
parser.add_argumane('--binary',action="store_true")
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
