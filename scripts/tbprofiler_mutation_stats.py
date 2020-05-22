import sys
import argparse
import csv
import os
from github import Github
import subprocess
import re
import statsmodels.api as sm
from collections import defaultdict, Counter
from tqdm import tqdm
import json
import numpy as np
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

def download_data():
    with open("/dev/null","w") as O:
        subprocess.call("wget http://tbdr.lshtm.ac.uk/VM/tbprofiler_results.tgz", shell=True, stderr=O, stdout=O)
        subprocess.call("tar -xvf tbprofiler_results.tgz", shell=True, stderr=O, stdout=O)

def get_codon_number(x):
    re_obj = re.search("p.[A-Za-z]+([0-9]+)[A-Za-z\*]+",x)
    return re_obj.group(1)

def main_identify_new_mutations(args):
    if not os.path.isfile("tbprofiler_results.tgz"):
        download_data()
    gene2locustag = {}
    drug2genes = defaultdict(set)
    conf = get_conf_dict(sys.base_prefix + "/share/tbprofiler/%s" % args.db)
    for l in open(conf["bed"]):
        row = l.rstrip().split()
        gene2locustag[row[4]] = row[3]
        gene2locustag[row[3]] = row[3]
        for d in row[5].split(","):
            drug2genes[d].add(row[3])

    mutations = set()
    for row in csv.DictReader(open(args.csv)):
        mutations.add((row["Drug"],gene2locustag[row["Gene"]],row["Mutation"]))

    multi_change_codons = defaultdict(list)
    for drug,gene,mutation in mutations:
        if "any_missense_codon_" in mutation:
            codon_num = mutation.replace("any_missense_codon_","")
            multi_change_codons[(gene,codon_num)].append(drug.lower())
    for drug,gene,mutation in mutations:
        if "any_missense_codon_" in mutation:
            codon_num = mutation.replace("any_missense_codon_","")
            multi_change_codons[(gene,codon_num)].append(drug.lower())

    meta = {}
    for row in csv.DictReader(open("tb.dst.csv")):
        meta[row["id"]] = row

    samples = [x.replace(".results.json","") for x in os.listdir("%s/" % args.dir) if x[-13:]==".results.json"]
    variants = defaultdict(lambda:defaultdict(list))
    mutation_types = defaultdict(dict)
    lineages = {}
    sys.stderr.write("Loading tb-profiler results\n")
    for s in tqdm(samples):
        tmp = json.load(open("%s/%s.results.json" % (args.dir,s)))
        lineages[s] = tmp["sublin"]
        for var in tmp["dr_variants"]:
            variants[var["locus_tag"]][var["change"]].append(s)
            if "large_deletion" in var["type"]:
                variants[var["locus_tag"]]["large_deletion"].append(s)
            elif  "frameshift" in var["type"]:
                variants[var["locus_tag"]]["frameshift"].append(s)
            if var["type"]=="missense":
                codon_num = get_codon_number(var["change"])
                if (var["locus_tag"],codon_num) in multi_change_codons and var["drug"] in multi_change_codons[(var["locus_tag"],codon_num)]:
                    variants[var["locus_tag"]]["any_missense_codon_"+codon_num].append(s)
            mutation_types[(var["locus_tag"],var["change"])] = var["type"]
        for var in tmp["other_variants"]:
            variants[var["locus_tag"]][var["change"]].append(s)
            if "large_deletion" in var["type"]:
                variants[var["locus_tag"]]["large_deletion"].append(s)
            elif  "frameshift" in var["type"]:
                variants[var["locus_tag"]]["frameshift"].append(s)
            if var["type"]=="missense":
                codon_num = get_codon_number(var["change"])
                if (var["locus_tag"],codon_num) in multi_change_codons and var["drug"] in multi_change_codons[(var["locus_tag"],codon_num)]:
                    variants[var["locus_tag"]]["any_missense_codon_"+codon_num].append(s)
            mutation_types[(var["locus_tag"],var["change"])] = var["type"]

    sys.stderr.write("Collected %s unique variants in %s genes\n" % (sum([len(variants[x]) for x in variants]),len(variants)))
    results = []

    sys.stderr.write("Calculating metrics for selected variants\n")
    for drug,gene,mutation in tqdm(mutations):
        if drug not in meta[samples[0]]: quit("%s not in meta" % drug)
        if gene not in variants: quit("%s not in genotype files" % gene)
        # sys.stderr.write("Calculating metrics for %s %s with %s\n" % (gene,mutation,drug))

        result = {"gene":gene,"drug":drug,"mutation":mutation}
        t = [
                [0.5,0.5],
                [0.5,0.5]
             ]
        for s in samples:
            if s not in meta: continue
            if meta[s][drug]=="1" and s in variants[gene][mutation]:         t[0][0]+=1
            if meta[s][drug]=="0" and s in variants[gene][mutation]:         t[0][1]+=1
            if meta[s][drug]=="1" and s not in variants[gene][mutation]:     t[1][0]+=1
            if meta[s][drug]=="0" and s not in variants[gene][mutation]:     t[1][1]+=1
        t2 = sm.stats.Table2x2(np.asarray(t))
        result["OR"] = t2.oddsratio
        result["OR_pval"] = t2.oddsratio_pvalue()
        result["RR"] = t2.riskratio
        result["RR_pval"] = t2.riskratio_pvalue()
        result["table"] = t
        result["variant_type"] = mutation_types[(gene,mutation)]
        result["num_samples"] = len(variants[gene][mutation])
        result["lineage"] = ", ".join(["%s:%s" % (x,y) for x,y in Counter([lineages[s] for s in variants[gene][mutation] if ";" not in lineages[s]]).most_common()])
        results.append(result)

    for i in tqdm(range(len(results))):
        if results[i]["OR"]>10 and results[i]["OR_pval"]<args.pval_cutoff and results[i]["RR"]>1 and results[i]["RR_pval"]<args.pval_cutoff:
            results[i]["confidence"] = "high"
        elif 5<results[i]["OR"]<=10 and results[i]["OR_pval"]<args.pval_cutoff and results[i]["RR"]>1 and results[i]["RR_pval"]<args.pval_cutoff:
            results[i]["confidence"] = "moderate"
        elif 1<results[i]["OR"]<=5 and results[i]["OR_pval"]<args.pval_cutoff and results[i]["RR"]>1 and results[i]["RR_pval"]<args.pval_cutoff:
            results[i]["confidence"] = "low"
        elif (results[i]["OR"]<=1 and results[i]["OR_pval"]<args.pval_cutoff) or (results[i]["RR"]<=1 and results[i]["RR_pval"]<args.pval_cutoff):
            results[i]["confidence"] = "no_association"
        else:
            results[i]["confidence"] = "indeterminate"

    fields = ["drug","gene","mutation","variant_type","table","num_samples","lineage","OR","OR_pval","RR","RR_pval","confidence"]
    with open(args.out,"w") as OUT:
        writer = csv.DictWriter(OUT, fields)
        writer.writeheader()
        for res in results:
            writer.writerow(res)

parser = argparse.ArgumentParser(description='Mutations stats',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--csv',help='CSV file containing Drug, Gene and Mutation fields',required=True)
parser.add_argument('--out',type=str,default="confidence.csv", help="Output file")
parser.add_argument('--dir',default="tbprofiler_results/",type=str,help='Firectory to look for tbprofiler results files')
parser.add_argument('--pval-cutoff',default=0.05,type=float,help='Pvalue cutoff to use for the corrected OR and RR p-vaule significance')
parser.add_argument('--db',default="tbdb",type=str,help='NGS Platform')
parser.set_defaults(func=main_identify_new_mutations)

args = parser.parse_args()
args.func(args)
