from collections import defaultdict
from pathogenprofiler import Vcf
import argparse
import csv
import logging
import re

def process_tb_profiler_args(args: argparse.Namespace) -> None:
    if args.snp_dist or args.update_phylo:
        args.call_whole_genome = True
    args.call_lineage = False if args.no_lineage else True
    if args.vcf and args.spoligotype:
        args.spoligotype = False
    if args.snp_dist or args.update_phylo:
        args.call_whole_genome = True

def get_vcf_samples(vcf_file):
    vcf = Vcf(vcf_file)
    return vcf.samples

def get_lt2drugs(bed_file):
    lt2drugs = {}
    for l in open(bed_file):
        row = l.strip().split()
        lt2drugs[row[3]] = row[5].split(",")
    return lt2drugs

def get_gene2drugs(bed_file):
    lt2drugs = {}
    for l in open(bed_file):
        row = l.strip().split()
        lt2drugs[row[4]] = row[5].split(",")
    return lt2drugs

def get_drugs2lt(bed_file):
    tmp = get_lt2drugs(bed_file)
    results = defaultdict(list)
    for gene in tmp:
        for drug in tmp[gene]:
            results[drug].append(gene)
    return dict(results)

def get_drugs2gene(bed_file):
    tmp = get_gene2drugs(bed_file)
    results = defaultdict(list)
    for gene in tmp:
        for drug in tmp[gene]:
            results[drug].append(gene)
    return dict(results)

def get_drug_list(bed_file):
    tmp = get_drugs2lt(bed_file)
    return set(tmp.keys())

def rv2genes(bed_file):
    #Chromosome      759310  763325  Rv0667  rpoB    rifampicin
    rv2gene = {}
    for l in open(bed_file):
        row = l.strip().split()
        rv2gene[row[3]] = row[4]
    return rv2gene

def genes2rv(bed_file):
    #Chromosome      759310  763325  Rv0667  rpoB    rifampicin
    rv2g = rv2genes(bed_file)
    gene2rv = {v:k for k,v in rv2g.items()}
    return gene2rv

def reformat_variant_csv_file(files: list, outfile: str) -> str:
    rows = []
    include_mutation = False
    for filename in files:
        for row in csv.DictReader(open(filename)):
            new_rows = {
                'Gene': row['Gene'],
            }
            if "Mutation" in row:
                new_rows['Mutation'] = row['Mutation']
                include_mutation = True
            info = {}
            for k,v in row.items():
                if k in ('Gene','Mutation'):
                    continue
                info[k] = v
            info_string = ";".join([f"{k.lower()}={v}" for k,v in info.items()])
            new_rows['Info'] = info_string
            rows.append(new_rows)


    with open(outfile, 'w') as csvfile:
        fieldnames = ['Gene','Mutation','Info'] if include_mutation else ['Gene','Info']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    
    return outfile

def check_db_version(db_version: str, tbprofiler_version: str) -> None:
    for d in db_version.split(","):
        r = re.search('([<>=]+)(.*)',d)
        if r==None:
            logging.error(f"Invalid version string: {d}")
            quit(1)

        d = f"{r.group(1)} '{r.group(2)}'"
        if eval(f"'{tbprofiler_version}' {d}")==False:
            if ">" in d:
                logging.error(f"Your version of tb-profiler ({tbprofiler_version}) is too old to use this version of the database. Please update tb-profiler to {db_version}")
            else:
                logging.error(f"Your version of tb-profiler ({tbprofiler_version}) is too new to use this version of the database. Please update the database to {db_version}")
            quit(1)