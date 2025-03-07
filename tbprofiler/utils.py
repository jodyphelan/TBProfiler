from collections import defaultdict
from pathogenprofiler import Vcf, get_db
import argparse
import csv
import logging
import re
from packaging.version import Version
import sys

def process_tb_profiler_args(args: argparse.Namespace) -> None:
    if args.snp_dist:
        args.call_whole_genome = True
    args.call_lineage = False if args.no_lineage else True
    if args.vcf and args.spoligotype:
        args.spoligotype = False

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

def check_db_version(db_current_version_str: str, compatible_schema_version_str: str):
    db_current_version = Version(db_current_version_str)
    compatible_schema_version = Version(compatible_schema_version_str)
    logging.debug(f"Database version: {db_current_version}")
    logging.debug(f"Compatible schema version: {compatible_schema_version}")
    if db_current_version.major != compatible_schema_version.major:
        logging.error(f"Latest database schema version {db_current_version_str} is not compatible with this version of tb-profiler (requires {compatible_schema_version.major}.x.x). Please make sure you are using the latest software and database versions.")
        quit(1)


def get_tier1_genes(dbname,dbdir=f'{sys.base_prefix}/share/tbprofiler'):
    conf = get_db(dbdir,dbname)
    id2name = rv2genes(conf['bed'])
    db = conf['json_db']
    tier1_genes = set()
    for gene in db:
        for var in db[gene]:
            for ann in db[gene][var]['annotations']:
                if ann['type']=='drug_resistance':
                    tier1_genes.add(id2name[gene])
    return tier1_genes

def get_default_db_dir():
    return f'{sys.base_prefix}/share/tbprofiler'