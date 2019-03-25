#! /usr/bin/env python
import sys
import argparse
import json
from collections import defaultdict
import re

indel_re = re.compile("([0-9]+)([A-Z]+)>([A-Z]+)")

def mutation_type(x):
	tmp = indel_re.search(x)
	if tmp:
		if len(tmp.group(2))!=len(tmp.group(3)):
			return "indel"
	return "snp"

def load_genes(bed):
	res = {}
	for l in open(bed):
		row = l.rstrip().split()
		res[row[3]] = row[4]
	return res

def load_library(library_file):
	lib = json.load(open(library_file))
	results = defaultdict(lambda: defaultdict(lambda: {"indels":[],"snps":[]}))
	for locus in lib:
		for var in lib[locus]:
			for drug in lib[locus][var]["drugs"]:
				if mutation_type(var)=="indel":
					results[drug][locus]["indels"].append(var)
				else:
					results[drug][locus]["snps"].append(var)
	return results

def main(args):
	library_file = "%s.dr.json" % args.prefix
	bed_file = "%s.bed" % args.prefix
	rv2gene = load_genes(bed_file)
	lib = load_library(library_file) #'Rv0682': {'indels': [], 'snps': ['86R>86P', '86R>86W', '9R>9H', '84G>84V', '43K>43R', '43K>43T', '51K>51N', '88K>88R', '88K>88Q', '88K>88M', '88K>88T', '40T>40I', '41T>41S', '52V>52G', '87V>87L', '93V>93M']
	print("Drug\tLocus_tag\tGene\tSNPs\tINDELs")
	drugs = [x.rstrip().lower() for x in open(args.drugs).readlines()] if args.drugs else list(lib.keys())
	if args.ngs:
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
	for drug in drugs:
		for locus in sorted(lib[drug]):
			if args.ngs:
				num_snps_ngs = sum([1 if v in variants[locus] for v in lib[drug][locus]["snps"]])
				num_indels_ngs = sum([1 if v in variants[locus] else 0 for v in lib[drug][locus]["indels"]])
				print("%s\t%s\t%s\t%s (%s)\t%s (%s)" % (drug,locus,rv2gene[locus],len(lib[drug][locus]["snps"]),num_snps_ngs,len(lib[drug][locus]["indels"]),num_indels_ngs))
			else:
				print("%s\t%s\t%s\t%s\t%s" % (drug,locus,rv2gene[locus],len(lib[drug][locus]["snps"]),len(lib[drug][locus]["indels"])))

def compare(args):
	library_file1 = "%s.dr.json" % args.prefix1
	library_file2 = "%s.dr.json" % args.prefix2
	lib1 = load_library(library_file1)
	lib2 = load_library(library_file2)
	new_drugs = {library_file1:set(),library_file2:set()}
	new_mutations = {"snps":{library_file1:set(),library_file2:set()},"indels":{library_file1:set(),library_file2:set()}}
	drugs = [x.rstrip().lower() for x in open(args.drugs).readlines()] if args.drugs else list(set(list(lib1.keys())).union(set(list(lib2.keys()))))
	lib1_loci = set()
	lib2_loci = set()
	for drug in drugs:
		if drug not in lib1 and drug not in lib2: continue
		loci = set()
		if drug in lib1:
			loci = loci.union(set(lib1[drug].keys()))
		else:
			new_drugs[library_file2].add(drug)
		if drug in lib2:
			loci = loci.union(set(lib2[drug].keys()))
		else:
			new_drugs[library_file1].add(drug)
		for locus in loci:
			#### snps ####
			lib1_set = set()
			if drug in lib1 and locus in lib1[drug]:
				lib1_set = set([json.dumps(x) for x in lib1[drug][locus]["snps"]])
				lib1_loci.add(locus)
			lib2_set = set()
			if drug in lib2 and locus in lib2[drug]:
				lib2_set = set([json.dumps(x) for x in lib2[drug][locus]["snps"]])
				lib2_loci.add(locus)
			for d in lib1_set.difference(lib2_set):
				new_mutations["snps"][library_file1].add((locus,d))
			for d in lib2_set.difference(lib1_set):
				new_mutations["snps"][library_file2].add((locus,d))

			lib1_set_indels = set()
			if drug in lib1 and locus in lib1[drug]:
				lib1_set_indels = set([json.dumps(x) for x in lib1[drug][locus]["indels"]])
			lib2_set_indels = set()
			if drug in lib2 and locus in lib2[drug]:
				lib2_set_indels = set([json.dumps(x) for x in lib2[drug][locus]["indels"]])
			for d in lib1_set_indels.difference(lib2_set_indels):
				new_mutations["indels"][library_file1].add((locus,d))
			for d in lib2_set_indels.difference(lib1_set_indels):
				new_mutations["indels"][library_file2].add((locus,d))
			print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (args.prefix1,args.prefix2,drug,locus,len(lib1_set.difference(lib2_set)),len(lib1_set.intersection(lib2_set)),len(lib2_set.difference(lib1_set)),len(lib1_set_indels.difference(lib2_set_indels)),len(lib1_set_indels.intersection(lib2_set_indels)),len(lib2_set_indels.difference(lib1_set_indels))))
			if locus=="Rv2043c":
				print(lib2_set_indels.difference(lib1_set_indels))
	print("drugs\t%s\t%s" % (len(new_drugs[library_file1]),len(new_drugs[library_file2])))
	print("loci\t%s\t%s" % (list(lib1_loci.difference(lib2_loci)),list(lib2_loci.difference(lib1_loci))))
	print("snps\t%s\t%s" % (len(new_mutations["snps"][library_file1]),len(new_mutations["snps"][library_file2])))
	print("indels\t%s\t%s" % (len(new_mutations["indels"][library_file1]),len(new_mutations["indels"][library_file2])))

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_sub = subparsers.add_parser('mutations', help='Run whole pipeline', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('prefix',help='NGS Platform')
parser_sub.add_argument('--drugs',default=None,help='NGS Platform')
parser_sub.add_argument('--ngs',default=None,help='NGS Platform')
parser_sub.set_defaults(func=main)

parser_sub = subparsers.add_parser('compare', help='Run whole pipeline', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('prefix1',help='NGS Platform')
parser_sub.add_argument('prefix2',help='NGS Platform')
parser_sub.add_argument('--drugs',default=None,help='NGS Platform')
parser_sub.set_defaults(func=compare)

args = parser.parse_args()
if vars(args)=={}:
	parser.print_help(sys.stderr)
else:
	args.func(args)
