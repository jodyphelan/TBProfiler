#! /usr/bin/env python
from __future__ import division
import pathogenprofiler as pp
import sys
import json
import csv
import argparse
from collections import defaultdict
from tqdm import tqdm

fluoroquinolones = ['ciprofloxacin','levofloxacin', 'moxifloxacin', 'ofloxacin','fluoroquinolones']
aminoglycosides = ['amikacin', 'capreomycin','kanamycin']
first_line = ["rifampicin","isoniazid","ethambutol","pyrazinamide"]


def write_itol(false_positives,false_negatives,drug):
	O = open("%s.fp_fn.itol.txt" % drug,"w")
	O.write("""DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL	%s
COLOR	#ff0000

LEGEND_TITLE	%s
LEGEND_SHAPES	1	1
LEGEND_COLORS	red	blue
LEGEND_LABELS	False_negative	False_positive

DATA
""" % (drug,drug))
	for x in false_positives:
		O.write("%s\tblue\n" % x)
	for x in false_negatives:
		O.write("%s\tred\n" % x)
	O.close()

def calculate(args):
	sample_file = args.samples
	dst_file = args.dst

	dst = json.load(open(dst_file))
	drug_loci = pp.load_bed(args.bed,[6],4) # {'Rv0668': ('rifampicin')}
	FAIL = open("samples_not_found.txt","w")
	samples = [x.rstrip() for x in open(sample_file).readlines()]
	ext = ".results.json"
	drugs = [d.lower() for d in dst[samples[0]].keys()]
	results = {d:{"tp":[],"tn":[],"fp":[],"fn":[]} for d in drugs+["flq","mdr","xdr","sus"]}
	counts = {d:{"tp":0,"tn":0,"fp":0,"fn":0} for d in drugs+["flq","mdr","xdr","sus"]}
	pre = args.dir if args.dir else ""
	for s in tqdm(samples):
		res_file = "%s/%s%s" % (pre,s,ext)
		if pp.nofile(res_file):
			pp.log("Warning: %s does not exist!" % res_file)
			FAIL.write("%s\n" % s)
			continue
		res = json.load(open(res_file))
		na_drugs = set()
		for locus in drug_loci:
			if res["missing_regions"][locus]>args.miss:
				for tmp in drug_loci[locus][0].split(","):
					na_drugs.add(tmp)
		resistant_drugs = [d["drug"].lower() for d in res["dr_variants"]]
		for d in drugs:
			if d in na_drugs:
				dst[s][d]="NA"

		for d in drugs:
			if dst[s][d]=="0" and d not in resistant_drugs:
				results[d]["tn"].append(s)
				counts[d]["tn"]+=1
			elif dst[s][d]=="0" and d in resistant_drugs:
				results[d]["fp"].append(s)
				counts[d]["fp"]+=1
			elif dst[s][d]=="1" and d not in resistant_drugs:
				results[d]["fn"].append(s)
				counts[d]["fn"]+=1
			elif dst[s][d]=="1" and d in resistant_drugs:
				results[d]["tp"].append(s)
				counts[d]["tp"]+=1

		#### Fluoroquinolones ####
		dst_flq = "0"
		dst_flq_NA = True

		for d in fluoroquinolones:
			if d not in dst[s]: continue
			if dst[s][d]!="NA": dst_flq_NA = False
			if dst[s][d]=="1": dst_flq = "1"

		dst_flq_list = [dst[s][d] for d in fluoroquinolones if d in dst[s]]
		if "1" in dst_flq_list and "0" in dst_flq_list:
			dst_flq = "NA"
		if dst_flq_NA: dst_flq = "NA"

		gst_flq = "0"
		for d in fluoroquinolones:
			if d in resistant_drugs: gst_flq = "1"

		if dst_flq=="1" and gst_flq=="1":
		  results["flq"]["tp"].append(s)
		  counts["flq"]["tp"]+=1
		if dst_flq=="0" and gst_flq=="1":
		  results["flq"]["fp"].append(s)
		  counts["flq"]["fp"]+=1
		if dst_flq=="1" and gst_flq=="0":
		  results["flq"]["fn"].append(s)
		  counts["flq"]["fn"]+=1
		if dst_flq=="0" and gst_flq=="0":
		  results["flq"]["tn"].append(s)
		  counts["flq"]["tn"]+=1

		#### MDR & XDR ####
		dst_mdr = "1" if dst[s]["rifampicin"]=="1" and dst[s]["isoniazid"]=="1" else "0"
		if dst[s]["rifampicin"]=="NA" or dst[s]["isoniazid"]=="NA": dst_mdr = "NA"
		flq = False
		flq_NA = True
		for d in fluoroquinolones:
			if d not in dst[s]: continue
			if dst[s][d]!="NA": flq_NA = False
			if dst[s][d]=="1": flq = True
		amg = False
		amg_NA = True
		for d in aminoglycosides:
			if d not in dst[s]: continue
			if dst[s][d]!="NA": amg_NA = False
			if dst[s][d]=="1": amg = True
		dst_xdr = "1" if dst_mdr=="1" and flq and amg else "0"

		if flq_NA or amg_NA: dst_xdr="NA"
		if dst_mdr=="NA": dst_xdr = "NA"

		#### Profiling results #####
		gst_mdr = "1" if "rifampicin" in resistant_drugs and "isoniazid" in resistant_drugs else "0"
		flq = False
		for d in fluoroquinolones:
			if d in resistant_drugs: flq = True
		amg = False
		for d in aminoglycosides:
			if d in resistant_drugs: amg = True
		gst_xdr = "1" if gst_mdr=="1" and flq and amg else "0"
		if dst_mdr=="1" and gst_mdr=="1":
			results["mdr"]["tp"].append(s)
			counts["mdr"]["tp"]+=1
		if dst_mdr=="0" and gst_mdr=="1":
			results["mdr"]["fp"].append(s)
			counts["mdr"]["fp"]+=1
		if dst_mdr=="1" and gst_mdr=="0":
			results["mdr"]["fn"].append(s)
			counts["mdr"]["fn"]+=1
		if dst_mdr=="0" and gst_mdr=="0":
			results["mdr"]["tn"].append(s)
			counts["mdr"]["tn"]+=1
		if dst_xdr=="1" and gst_xdr=="1":
		  results["xdr"]["tp"].append(s)
		  counts["xdr"]["tp"]+=1
		if dst_xdr=="0" and gst_xdr=="1":
		  results["xdr"]["fp"].append(s)
		  counts["xdr"]["fp"]+=1
		if dst_xdr=="1" and gst_xdr=="0":
		  results["xdr"]["fn"].append(s)
		  counts["xdr"]["fn"]+=1
		if dst_xdr=="0" and gst_xdr=="0":
		  results["xdr"]["tn"].append(s)
		  counts["xdr"]["tn"]+=1
		### susceptibility
		if "NA" not in [dst[s][d] for d in first_line]:
			dst_sus = "1" if "1" not in [dst[s][d] for d in drugs] else "0"
			gst_sus = "1" if all([x not in resistant_drugs for x in first_line]) else "0"
			if dst_sus=="1" and gst_sus=="1":
			  results["sus"]["tp"].append(s)
			  counts["sus"]["tp"]+=1
			if dst_sus=="0" and gst_sus=="1":
			  results["sus"]["fp"].append(s)
			  counts["sus"]["fp"]+=1
			if dst_sus=="1" and gst_sus=="0":
			  results["sus"]["fn"].append(s)
			  counts["sus"]["fn"]+=1
			if dst_sus=="0" and gst_sus=="0":
			  results["sus"]["tn"].append(s)
			  counts["sus"]["tn"]+=1
	json.dump(results,open("results.json","w"))
	json.dump(counts,open("counts.json","w"))
	counts = json.load(open("counts.json"))
	drugs = [x.rstrip().lower() for x in open(args.drugs).readlines()] if args.drugs else list(counts.keys())
	print("Drug\tNum\tSusceptible\tResistant\tSensitivity\tSpecificity")
	for d in drugs:
		if d not in counts: continue
		if counts[d]["tp"]+counts[d]["fn"]==0 or counts[d]["tn"]+counts[d]["fp"]==0: continue
		sensitivity = counts[d]["tp"]/(counts[d]["tp"]+counts[d]["fn"])
		specificity = counts[d]["tn"]/(counts[d]["tn"]+counts[d]["fp"])
		total = counts[d]["tp"]+counts[d]["fp"]+counts[d]["tn"]+counts[d]["fn"]
		suc = counts[d]["tn"]+counts[d]["fp"]
		res = counts[d]["tp"]+counts[d]["fn"]
		print("%s\t%s\t%s\t%s\t%s\t%s" % (d.capitalize(),total,suc,res,sensitivity,specificity))

def print_numbers(args):
	counts = json.load(open("counts.json"))
	drugs = [x.rstrip().lower() for x in open(args.drugs).readlines()] if args.drugs else list(counts.keys())
	print("Drug\tNum\tSusceptible\tResistant\tSensitivity\tSpecificity\tAccuracy\tPPV\tNPV")
	for d in drugs:
		if d not in counts: continue
		if counts[d]["tp"]+counts[d]["fn"]==0 or counts[d]["tn"]+counts[d]["fp"]==0: continue
		sensitivity = counts[d]["tp"]/(counts[d]["tp"]+counts[d]["fn"])
		specificity = counts[d]["tn"]/(counts[d]["tn"]+counts[d]["fp"])
		accuracy = (counts[d]["tp"]+counts[d]["tn"])/(counts[d]["tn"]+counts[d]["tp"]+counts[d]["fn"]+counts[d]["fp"])
		ppv = counts[d]["tp"]/(counts[d]["tp"]+counts[d]["fp"]) if (counts[d]["tp"]+counts[d]["fp"]) > 0 else 0.0
		npv = counts[d]["tn"]/(counts[d]["tn"]+counts[d]["fn"]) if (counts[d]["tn"]+counts[d]["fn"]) > 0 else 0.0
		total = counts[d]["tp"]+counts[d]["fp"]+counts[d]["tn"]+counts[d]["fn"]
		suc = counts[d]["tn"]+counts[d]["fp"]
		res = counts[d]["tp"]+counts[d]["fn"]
		print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (d.capitalize(),total,suc,res,sensitivity,specificity,accuracy,ppv,npv))

def analyse(args):
	drug_loci = pp.load_bed(args.bed,[6],4)
	for d in drug_loci:
		drug_loci[d] = [x.lower() for x in drug_loci[d][0].split(";")]

	drug =	args.drug
	results = json.load(open("results.json"))
	fn = defaultdict(int)
	fp = defaultdict(int)
	print("Number true positives: %s" % len(results[drug]["tp"]))
	print("Number false negatives: %s" % len(results[drug]["fn"]))
	print("Number true negatives: %s" % len(results[drug]["tn"]))
	print("Number false positives: %s" % len(results[drug]["fp"]))
	for s in results[drug]["fn"]:
		print("%s\tFalse_negative" % s)
		tmp = json.load(open("%s/%s.results.json"%(args.dir,s)))
		for var in tmp["other_variants"]:
			if drug not in drug_loci[var["locus_tag"]]:
				continue
			fn[(var["gene"],var["locus_tag"],var["change"])]+=1
	for s in results[drug]["fp"]:
		print("%s\tFalse_positive" % s)
		tmp = json.load(open("%s/%s.results.json"%(args.dir,s)))
		for var in tmp["dr_variants"]:
			#if drug not in drug_loci[var["locus_tag"]]:
			if drug != var["drug"].lower():
				continue

			fp[(var["gene"],var["locus_tag"],var["change"])]+=1
	for key in sorted(fn,key=lambda x:fn[x]):
		gene,rv,var = key
		print("False_Negative\t%s\t%s\t%s" % (gene,var,fn[key]))
	for key in sorted(fp,key=lambda x:fp[x]):
		gene,rv,var = key
#		tmp = odds_ratio[rv][var]
#		print "False_Positive\t%s\t%s\t%s\t%s\t%s\t%s" % (gene,var,fp[key],tmp[4],tmp[5],tmp[6])
		print("False_Positive\t%s\t%s\t%s" % (gene,var,fp[key]))
	if args.itol:
		write_itol(results[drug]["fp"],results[drug]["fn"],args.drug)

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_sub = subparsers.add_parser('calculate', help='Run whole pipeline', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('samples',help='NGS Platform')
parser_sub.add_argument('dst',help='NGS Platform')
parser_sub.add_argument('bed',help='NGS Platform')
parser_sub.add_argument('--dir','-d',default="results/",type=str,help='NGS Platform')
parser_sub.add_argument('--miss','-m',default=0.1,type=float,help='Fraction of gene missing to call gDST as missing')
parser_sub.add_argument('--drugs',default=None,type=str,help='NGS Platform')
parser_sub.set_defaults(func=calculate)

parser_sub = subparsers.add_parser('print', help='Run whole pipeline', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--drugs','-d',default=None,type=str,help='NGS Platform')
parser_sub.set_defaults(func=print_numbers)

parser_sub = subparsers.add_parser('analyse', help='Run whole pipeline', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('drug',help='NGS Platform')
parser_sub.add_argument('bed',help='NGS Platform')
parser_sub.add_argument('--dir','-d',default="results/",type=str,help='NGS Platform')
parser_sub.add_argument('--itol',action="store_true",help='NGS Platform')
parser_sub.set_defaults(func=analyse)

args = parser.parse_args()
if vars(args)=={}:
	parser.print_help(sys.stderr)
else:
	args.func(args)
