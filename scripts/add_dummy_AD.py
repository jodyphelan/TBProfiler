#! /usr/bin/env python3
import sys
import re
import argparse

def main(args):
	if args.ref:
		chroms = []
		for l in open(args.ref+".fai"):
			chroms.append(l.strip().split())
	for l in sys.stdin:
		row = l.rstrip().split()
		if l[0]=="#":
			# temp fix for issue with pilon vcf header
			if "ID=DP,Number=1,Type=String" in l:
				l = re.sub(r'ID=DP,Number=1,Type=String', r'ID=DP,Number=1,Type=Integer',l)
			if "source=lofreq call" in l:
				r = re.search("(\S+).bam",l)
				sample_name = args.sample_name if args.sample_name else r.group(1) 
				sys.stdout.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
				sys.stdout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
				for chrom in chroms:
					sys.stdout.write(f'##contig=<ID={chrom[0]},length={chrom[1]}>\n')
			if "#CHROM" in l:
				sys.stdout.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)">\n')
				if len(row)==8:
					l = l.strip() + f"\tFORMAT\t{sample_name}\n"
			sys.stdout.write(l)
		else:
			if len(row)<9:
				row.append("GT")
				row.append("1/1")
			if "AD" not in row[8]:
				row[8]+=":AD"
				if "DP4=" in row[7]:
					r = re.search("DP4=([0-9,]+)",row[7])
					if r:
						dp4 = [int(x) for x in r.group(1).split(",")]
						ref = dp4[0]+dp4[1]
						alt = dp4[2]+dp4[3]
						row[9]+=f":{ref},{alt}"
				elif "AF=" in row[7]:
					r = re.search("AF=([0-9\.]+)",row[7])
					if r:
						alt = int(100*float(r.group(1)))
						ref = int(100-alt)
						row[9]+=f":{ref},{alt}"
				
				else:
					alt = 100
					ref = 100
					row[9]+=f":{ref},{alt}"
			if "DP" not in row[8]:
				row[8]+=":DP"
				row[9]+=f":{ref+alt}"
		
			sys.stdout.write("%s\n" % "\t".join(row))

parser = argparse.ArgumentParser(description='add required annotations',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--ref',type=str,help='')
parser.add_argument('--sample-name',type=str,help='')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)