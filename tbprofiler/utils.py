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



import pysam
from scipy.stats import ranksums



rc = {"A": "T", "T": "A", "C": "G", "G": "C"}

def read_pos_rank_sum(bam: str, chrom: str, pos: int, ref: str, alt: str):
    bam = pysam.AlignmentFile(bam, "rb")
    
    ref_positions = []
    alt_positions = []
    
    for read in bam.fetch(chrom, pos, pos + 1):
        if read.is_unmapped or read.is_duplicate or read.is_secondary:
            continue
        
        # Get the read position that aligns to the genomic position
        read_pos = read.get_reference_positions(full_length=True)
        try:
            read_idx = read_pos.index(pos-1)
        except ValueError:
            continue  # This read doesn't actually cover the position (e.g., deletion)

        base = read.query_sequence[read_idx]
        
        # Store read-relative position (0 to read length)
        rel_pos = read_idx
        

        if base.upper() == ref.upper():
            ref_positions.append(rel_pos)
        elif base.upper() == alt.upper():
            alt_positions.append(rel_pos)
        else:
            continue  # Ignore mismatches that don't match ref or alt
    

    if len(ref_positions) == 0 or len(alt_positions) == 0:
        print("Insufficient data: need both ref and alt reads.")
        return None

    # Perform Wilcoxon rank-sum test
    stat, p_value = ranksums(alt_positions, ref_positions)
    
    return {
        "statistic": stat,
        "p_value": p_value,
        "alt_positions": alt_positions,
        "ref_positions": ref_positions
    }

def alt_near_read_ends_pct(bam: str, chrom: str, pos: int, alt: str, margin: int = 10):
    bam = pysam.AlignmentFile(bam, "rb")
    
    alt_read_positions = []
    near_end_count = 0
    
    for read in bam.fetch(chrom, pos, pos + 1):
        if read.is_unmapped or read.is_duplicate or read.is_secondary:
            continue
        
        # Get aligned reference-to-read positions
        read_pos = read.get_reference_positions(full_length=True)
        try:
            read_idx = read_pos.index(pos-1)
        except ValueError:
            continue  # skip reads that do not align to the position (e.g. deletions)

        # Skip if out of read bounds
        if read_idx >= len(read.query_sequence):
            continue

        base = read.query_sequence[read_idx]

        
        if base.upper() == alt.upper():
            alt_read_positions.append(read_idx)
            read_len = read.query_length
            # Check if variant is within 'margin' bases from either end
            if read_idx < margin or (read_len - read_idx - 1) < margin:
                near_end_count += 1

    total_alt_reads = len(alt_read_positions)
    if total_alt_reads == 0:
        return {
            "percent_near_ends": None,
            "total_alt_reads": 0,
            "near_end_reads": 0
        }

    percent = 100.0 * near_end_count / total_alt_reads
    return {
        "percent_near_ends": percent,
        "total_alt_reads": total_alt_reads,
        "near_end_reads": near_end_count
    }

import math

def calculate_epp(bam: str, chrom: str, pos: int, alt: str, margin: int = 10):
    bam = pysam.AlignmentFile(bam, "rb")
    EL = 0
    ER = 0

    for read in bam.fetch(chrom, pos, pos + 1):
        # check read name
        if read.is_unmapped or read.is_duplicate or read.is_secondary:
            continue

        ref_pos = read.get_reference_positions(full_length=True)
        try:
            read_idx = ref_pos.index(pos-1)
        except ValueError:
            continue  # read doesn't align over this position

        if read_idx >= len(read.query_sequence):
            continue

        
        

        base = read.query_sequence[read_idx]

        if base.upper() != alt.upper():
            continue

        read_len = read.query_length
        if read_idx < margin:
            EL += 1
        elif (read_len - read_idx - 1) < margin:
            ER += 1


    n = EL + ER
    if n == 0:
        return {
            "EPP": None,
            "EL": EL,
            "ER": ER
        }

    delta = abs((EL / n) - 0.5)
    p = 2 * math.exp(-2 * n * delta ** 2)
    p = min(max(p, 1e-300), 1.0)  # Clamp for numerical safety

    epp = -10 * math.log10(p)
    return {
        "EPP": epp,
        "EL": EL,
        "ER": ER
    }




