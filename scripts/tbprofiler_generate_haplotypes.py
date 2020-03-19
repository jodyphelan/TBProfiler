#! /usr/bin/env python
import json
from collections import defaultdict
import argparse
import os
from tqdm import tqdm
import sys
import pathogenprofiler as pp


try:
    sys.base_prefix
except:
    sys.base_prefix = getattr(sys, 'base_prefix', getattr(sys, 'real_prefix', sys.prefix))


def get_conf_dict(library_prefix):
    files = {"gff":".gff","ref":".fasta","ann":".ann.txt","barcode":".barcode.bed","bed":".bed","json_db":".dr.json","version":".version.json"}
    conf = {}
    for key in files:
        sys.stderr.write("Using %s file: %s\n" % (key,library_prefix+files[key]))
        conf[key] = pp.filecheck(library_prefix+files[key])
    return conf



def main(args):
	change_field = "change" if args.variant_format=="hgvs" else "_internal_change"
	conf = conf = get_conf_dict(sys.base_prefix + "/share/tbprofiler/%s" % args.db)
	drug2genes = defaultdict(set)
	gene2drugs = defaultdict(set)
	for l in open(conf["bed"]):
		row = l.rstrip().split()
		for d in row[5].split(","):
			drug2genes[d].add(row[3])
			gene2drugs[row[3]].add(d)
	if args.samples:
		samples = [x.rstrip() for x in open(args.samples).readlines()]
	else:
		samples = [x.replace(".results.json","") for x in os.listdir(args.dir) if x[-13:]==".results.json"]
	variants = defaultdict(lambda:defaultdict(list))
	if args.non_dr:
		print("sample,%s" % (",".join(["%s,%s" % ("dr_mutations_%s" % x,"other_mutations_%s" % x)  for x in sorted(drug2genes)])))
	else:
		print("sample,%s" % (",".join(["dr_mutations_%s" % x for x in sorted(drug2genes)])))
	for s in tqdm(samples):
		if pp.nofile("%s/%s.results.json" % (args.dir,s)):
			if args.non_dr:
				print("%s,%s" % (s,",".join(["NA,NA" for _ in drug2genes])))
			else:
				print("%s,%s" % (s,",".join(["NA" for _ in drug2genes])))
			continue
		tmp = json.load(open("%s/%s.results.json" % (args.dir,s)))
		sample_dr_mutations = defaultdict(list)
		tmp_store = set([json.dumps(x) for x in tmp["dr_variants"]])	# This step is only to remove duplicate
		tmp["dr_variants"] = []											# mutations introduced by a bug that is
		tmp["dr_variants"] = [json.loads(x) for x in tmp_store]			# now fixed
		for var in tmp["dr_variants"]:
			sample_dr_mutations[var["drug"]].append(var["gene"]+"_"+var[change_field])
		if args.non_dr:
			sample_other_mutations = defaultdict(list)
			tmp_store = set([json.dumps(x) for x in tmp["other_variants"]])	# This step is only to remove duplicate
			tmp["other_variants"] = []										# mutations introduced by a bug that is
			tmp["other_variants"] = [json.loads(x) for x in tmp_store]		# now fixed
			for var in tmp["other_variants"]:
				for d in gene2drugs[var["locus_tag"]]:
					sample_other_mutations[d].append(var["gene"]+"_"+var[change_field])
			print("%s,%s" % (s,",".join(["%s,%s" % ("; ".join(sample_dr_mutations[d]) if d in sample_dr_mutations else "WT","; ".join(sample_other_mutations[d]) if d in sample_other_mutations else "WT",) for d in sorted(drug2genes)])))
		else:
			print("%s,%s" % (s,",".join(["; ".join(sample_dr_mutations[d]) if d in sample_dr_mutations else "WT" for d in sorted(drug2genes)])))

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--samples',type=str,help='NGS Platform')
parser.add_argument('--dir',default="results/",type=str,help='NGS Platform')
parser.add_argument('--db',default="tbdb",type=str,help='NGS Platform')
parser.add_argument('--variant-format',default="hgvs",choices=["hgvs","bcftools"],type=str,help='NGS Platform')
parser.add_argument('--non-dr',action="store_true",help='NGS Platform')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
