#! /usr/bin/env python
import json
import sys
from collections import defaultdict
import re

drugs = [
	'isoniazid', 'rifampicin', 'ethambutol', 'pyrazinamide', 'streptomycin',
	'fluoroquinolones', 'amikacin', 'capreomycin', 'kanamycin',
	 'cycloserine',  'ethionamide', 'clofazimine', 'para-aminosalicylic_acid',
	 'delamanid', 'bedaquiline', 'linezolid', ]

variants = defaultdict(lambda:defaultdict(int))
data = json.load(open(sys.argv[1]))
for s in data:
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
