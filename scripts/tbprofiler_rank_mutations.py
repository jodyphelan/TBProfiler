#! /usr/bin/env python
import json
from collections import defaultdict
import argparse
import os
from tqdm import tqdm
import sys
import pathogenprofiler as pp
from collections import Counter
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
   # conf = get_conf_dict(sys.base_prefix + "/share/tbprofiler/%s" % args.db)
    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(".results.json", "") for x in os.listdir(args.dir) if x[-13:] == ".results.json"]
    variants = []
    for s in samples:
        data_to_add = set()
        tmp = json.load(open(f"{args.dir}/{s}.results.json"))
        for var in tmp["dr_variants"]:
            nt = var["nucleotide_change"] if "nucleotide_change" in var else var["change"]
            data_to_add.add((nt,var["change"],var["gene"]))
        variants += list(data_to_add)
    count = Counter(variants)
    for c in count.most_common(len(count)):
        nt,aa,gene = c[0]
        print("%s\t%s\t%s\t%s" % (gene,nt,aa,c[1]))

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--samples',type=str,help='NGS Platform')
parser.add_argument('--dir',default="results/",type=str,help='NGS Platform')
parser.add_argument('--db',default="tbdb",type=str,help='NGS Platform')
parser.add_argument('--variant-format',default="hgvs",choices=["hgvs","bcftools"],type=str,help='NGS Platform')
parser.add_argument('--non-dr',action="store_true",help='NGS Platform')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
