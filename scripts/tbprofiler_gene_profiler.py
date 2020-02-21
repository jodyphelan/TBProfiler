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


def main_profile(args):
    if args.external_db:
        conf = get_conf_dict(args.external_db)
    else:
        conf = get_conf_dict(sys.base_prefix+"/share/tbprofiler/%s" % args.db)
    if not args.prefix:
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
    pp.run_cmd("rm %(prefix)s.targets.vcf.gz* %(prefix)s.targets.csq.vcf.gz*" % vars(args))


def main_collate(args):
    samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]
    variants = defaultdict(lambda:defaultdict(set))
    for s in samples:
        j = json.load(open("%s/%s%s" % (args.dir,s,args.suffix)))
        print("%s\t%s" % (s,";".join([var["change"] for var in j["variants"]])))
        for var in j["variants"]:
            variants[var["locus_tag"]][var["change"]].add(s)

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")



parser_sub = subparsers.add_parser('profile', help='(BETA) Run profiling pipeline on Fasta file. Warning: this assumes that this is a good quality assembly which coveres all drug resistance loci', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--bam',type=str,help='BAM file',required=True)
parser_sub.add_argument('--prefix',type=str,help='Sample prefix for all results generated')
parser_sub.add_argument('--suffix',type=str,default=".results.json",help='Output file suffix')
parser_sub.add_argument('--platform','-m',choices=["illumina","nanopore"],default="illumina",help='NGS Platform used to generate data')
parser_sub.add_argument('--caller',default="bcftools", choices=["bcftools","gatk"],help="Variant calling tool to use.",type=str)
parser_sub.add_argument('--db',default='tbdb',help='Mutation panel name')
parser_sub.add_argument('--external_db',type=str,help='Path to db files prefix (overrides "--db" parameter)')
parser_sub.add_argument('--threads','-t',default=1,help='Threads to use',type=int)
parser_sub.set_defaults(func=main_profile)

parser_sub = subparsers.add_parser('collate', help='Output program version and exit', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--dir',type=str,default=".",help='Result directory')
parser_sub.add_argument('--suffix',type=str,default=".results.json",help='Output file suffix')
parser_sub.set_defaults(func=main_collate)

args = parser.parse_args()
if vars(args)=={}:
    parser.print_help(sys.stderr)
else:
    args.func(args)
