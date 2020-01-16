#! /usr/bin/env python
import os
import sys
import argparse
import json
from collections import defaultdict
from tqdm import tqdm
import pathogenprofiler as pp


def get_conf_dict(library_prefix):
    files = {"gff":".gff","ref":".fasta","ann":".ann.txt","barcode":".barcode.bed","bed":".bed","json_db":".dr.json","version":".version.json"}
    conf = {}
    for key in files:
        sys.stderr.write("Using %s file: %s\n" % (key,library_prefix+files[key]))
        conf[key] = pp.filecheck(library_prefix+files[key])
    return conf


def main(args):
    if args.external_db:
        conf = get_conf_dict(args.external_db)
    else:
        conf = get_conf_dict(sys.base_prefix+"/share/tbprofiler/%s" % args.db)
    args.prefix  = args.bam.split("/")[-1].replace(".bam","").replace(".cram","")
    bam_obj = pp.bam(args.bam, args.prefix, platform=args.platform)
    vcf_obj = bam_obj.call_variants(conf["ref"], caller=args.caller, bed_file=conf["bed"], threads=args.threads)
    csq_vcf_obj = vcf_obj.csq(conf["ref"],conf["gff"])
    csq = csq_vcf_obj.load_csq(ann_file=conf["ann"])
    results = {"variants":[]}
    for sample in csq:
        results["variants"]  = csq[sample]
    outfile = "%s%s" % (args.prefix,args.suffix)
    json.dump(results,open(outfile,"w"))

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bam',type=str,help='NGS Platform',required=True)
parser.add_argument('--suffix',type=str,default=".results.json",help='NGS Platform')
parser.add_argument('--platform','-m',choices=["illumina","nanopore"],default="illumina",help='NGS Platform used to generate data')
parser.add_argument('--caller',default="bcftools", choices=["bcftools","gatk"],help="Variant calling tool to use.",type=str)
parser.add_argument('--db',default='tbdb',help='Mutation panel name')
parser.add_argument('--external_db',type=str,help='Path to db files prefix (overrides "--db" parameter)')
parser.add_argument('--threads','-t',default=1,help='Threads to use',type=int)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
