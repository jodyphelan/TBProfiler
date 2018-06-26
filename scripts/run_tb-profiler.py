#! /usr/bin/env python
import sys
import csv
import argparse
import os.path

def filecheck(filename):
	"""
	Check if file is there and quit if it isn't
	"""
	if not os.path.isfile(filename):
		print("Can't find %s" % filename)
		exit(1)
	else:
		return filename

def main(args):
	out_script = "tb-profiler.run.sh"
	O = open(out_script,"w")
	samples = []
	for row in csv.DictReader(open(args.sample_file)):
		paired = True
		params = {}
		params["r1"] = "%s/%s" % (args.fastq_dir,row["R2"])
		filecheck(params["r1"])
		params["r2"] = "%s/%s" % (args.fastq_dir,row["R1"])
		if params["r2"]=="NA":
			paired = False
		else:
			filecheck(params["r2"])
		params["prefix"] = row["ID"]
		params["threads"] = args.threads
		params["mapper"] = args.mapper
		params["platform"] = args.platform
		params["caller"] = args.caller
		params["verbosity"] = args.verbosity
		params["tb_profiler"] = args.tb_profiler_dir if args.tb_profiler_dir else "tb-profiler"
		params["read_str"] = "-1 %(r1)s -2 %(r2)s" % params  if paired else "-1 %(r1)s" % params
		O.write("%(tb_profiler)s --prefix %(prefix)s  %(read_str)s --platform %(platform)s --mapper %(mapper)s --threads %(threads)s --caller %(caller)s --verbose %(verbosity)s") % params
	O.close()


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('sample_file', help='CSV file containig the following columns: ID, R1 and R2 (ID of the sample and the two read names)')
parser.add_argument('--platform','-m',choices=["illumina","minION"],default="illumina", help='Sequencing platform')
parser.add_argument('--fastq_dir','-f',default=".",type=str, help='Directory containing fastqs')
parser.add_argument('--threads',"-t",type=int,default=1, help='Number of available threads')
parser.add_argument('--mapper',type=str,choices=["bwa","minimap2","bowtie2"],default="bwa", help='Mapping tool to use')
parser.add_argument('--caller',type=str,choices=["lofreq","bcftools"],default="lofreq", help='Variant caller to use')
parser.add_argument('--tb_profiler_dir',type=str,default=None, help='Full path to the tb-profiler script (use if tb-profiler is not in the path)')
parser.add_argument('--verbosity',type=int,choices=[0,1,2],default=0, help='Variant caller to use')

parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
