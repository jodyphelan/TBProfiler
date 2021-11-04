import pathogenprofiler as pp
import sys
import argparse


def main(args):
    generator = pp.cmd_out(f"bcftools view {args.vcf}") if args.vcf else sys.stdin
    convert = dict(zip(args.source,args.target))
    for l in generator:    
        if l[0]=="#":
            sys.stdout.write(l)
        else:
            row = l.strip().split()
            row[0] = convert[row[0]]
            sys.stdout.write("\t".join(row)+"\n")

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',type=str,help='')
parser.add_argument('--source',nargs="+",type=str,help='')
parser.add_argument('--target',nargs="+",type=str,help='')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
