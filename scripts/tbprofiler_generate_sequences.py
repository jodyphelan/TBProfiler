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
    conf = get_conf_dict(sys.base_prefix + "/share/tbprofiler/%s" % args.db)
    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(".targets.csq.vcf.gz", "") for x in os.listdir(args.dir) if x[-19:] == ".targets.csq.vcf.gz"]
    sample_fastas = defaultdict(list)
    params = {"tmp_locations": pp.get_random_file(), "tmp_mappings": pp.get_random_file(), "ref":conf["ref"]}
    pp.run_cmd("awk '{print $1\":\"$2\"-\"$3\"\\t\"$5}' %s > %s" % (conf["bed"],params["tmp_mappings"]))
    pp.run_cmd("cut -f1 %s > %s" % (params["tmp_mappings"],params["tmp_locations"]))
    FILES = {}
    for l in open(params["tmp_mappings"]):
        row = l.rstrip().split()
        FILES[row[0]] = open("%s.fasta" % row[1],"w")
    for s in samples:
        params["vcf"] = "%s/%s.targets.csq.vcf.gz" % (args.dir, s)
        params["sample_fa"] = "%s.targets.fa" % (s)
        pp.run_cmd("samtools faidx -r %(tmp_locations)s %(ref)s | bcftools consensus -H A %(vcf)s > %(sample_fa)s" % params)
        fa_dict = pp.fasta(params["sample_fa"]).fa_dict
        for locus in fa_dict:
            FILES[locus].write(">%s\n%s\n" % (s,fa_dict[locus]))
    pp.rm_files([params["tmp_locations"], params["tmp_mappings"]])

parser = argparse.ArgumentParser(
    description='TBProfiler pipeline',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
parser.add_argument('--samples', type=str, help='NGS Platform')
parser.add_argument('--dir', default = "vcf/", type = str, help = 'NGS Platform')
parser.add_argument('--db', default="tbdb", type=str, help='NGS Platform')
parser.add_argument('--variant-format', default="hgvs", choices=["hgvs", "bcftools"], type=str,help='NGS Platform')

parser.add_argument('--non-dr', action="store_true", help='NGS Platform')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
