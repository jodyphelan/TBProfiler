#! /home/jody/software/anaconda2/bin/python
import json
import sys
import subprocess
from collections import defaultdict
import argparse
import os

################### Literals #####################
scriptDir = os.path.dirname(os.path.realpath(__file__))
tabix = scriptDir+"/bin/tabix"
bgzip = scriptDir+"/bin/bgzip"
snap = scriptDir+"/bin/snap-aligner"
htsbox = scriptDir+"/bin/htsbox"
sambamba = scriptDir+"/bin/sambamba"
dr_bed_file = scriptDir+"/dr.bed"
ref_dir = scriptDir+"/ref"
ref_file = ref_dir+"/MTB-h37rv_asm19595v2-eg18.fa"
ref_annotation = ref_dir+"/MTB-h37rv_asm19595v2-eg18.tab.ann.gz"
################### Functions ####################

def check_files(r1,r2):
	programs = {"tabix":tabix,"bgzip":bgzip,"sambamba":sambamba,"snap":snap,"htsbox":htsbox}
	for p in programs:
		if not os.path.isfile(programs[p]):
			print "Can't find %s" % p
			quit()
	if not r1:
		print "Please provide at least one read file"
		quit()

	if not os.path.isfile(r1):
		print "Can't find %s" % r1
		quit()
	if r2:
		if not os.path.isfile(r2):
			print "Can't find %s" % r2
			quit()


def index_file(infile):
	import os
	import subprocess
	tbi = infile+".tbi"
	if not os.path.isfile(tbi):
		subprocess.call("%s -b 2 -e 2 -s 1 %s" % (tabix,infile),shell=True)




def load_ann(bed_name):
	ann_dict = defaultdict(dict)
	annCMD = "%s %s -R %s" % (tabix,ref_annotation,bed_name)
	annPIPE = subprocess.Popen(annCMD,shell=True,stdout=subprocess.PIPE)
	for l in annPIPE.stdout:
		arr = l.rstrip().split()
		nuc_obj = {}
		nuc_obj[arr[2]] = {"codon":arr[6],"aa":arr[11]}
		nuc_obj[arr[3]] = {"codon":arr[8],"aa":arr[12]}
		nuc_obj[arr[4]] = {"codon":arr[9],"aa":arr[13]}
		nuc_obj[arr[5]] = {"codon":arr[10],"aa":arr[14]}
		ann_dict[arr[0]][arr[1]] = {"ref_nt":arr[2],"ref_codon":arr[6],"ref_aa":arr[11],"chr":arr[0],"pos":arr[1],"ann":nuc_obj,"rv":arr[15],"gene":arr[16],"gene_syn":arr[17],"ncr":arr[18],"start":arr[19],"end":arr[20],"strand":arr[21],"drug":arr[22],"ppe":arr[23],"codon_num":arr[24],"gene_nt":arr[25],"operon":arr[26]}
	return ann_dict

def load_bed(bed_name):
	temp_bed = defaultdict(lambda : defaultdict(set))
	for l in open(bed_name):
		arr = l.rstrip().split()
		drugs = arr[3].split(";")
		for d in drugs:
			temp_bed[arr[0]][arr[1]].add(d)
	return temp_bed

def doMapping(r1,r2,threads,prefix):
	if r2:
		cmd_mapping = "%s paired %s -compressedFastq %s %s -t %s -o -bam - | %s sort /dev/stdin -o %s.bam -t %s" % (snap,ref_dir,r1,r2,threads,sambamba,prefix,threads)
	else:
		cmd_mapping = "%s single %s -compressedFastq %s -t %s -o -bam - | %s sort /dev/stdin -o %s.bam -t %s" % (snap,ref_dir,r1,threads,sambamba,prefix,threads)
	subprocess.call(cmd_mapping,shell=True)

def doPileup(prefix):
	cmd_pileup = "%s pileup -f %s %s.bam | %s -c > %s.pileup.gz" % (htsbox,ref_file,prefix,bgzip,prefix)
	subprocess.call(cmd_pileup,shell=True)

def getCalls(prefix,min_cov,read_frac):
	pileup_file = prefix+".pileup.gz"
	base_calls = defaultdict(lambda : defaultdict(dict))
	index_file(pileup_file)
	for l in open(dr_bed_file):
		arr = l.rstrip().split()
		base_calls[arr[0]][arr[1]] = "NA"
	cmd_tabix = "%s %s -T %s" % (tabix,pileup_file,dr_bed_file)
	tabix_out = subprocess.Popen(cmd_tabix,shell=True,stdout=subprocess.PIPE)
	for l in tabix_out.stdout:
		chrom,pos,ref,temp_alts,temp_meta =  l.rstrip().split()
		cov = temp_meta.split(":")[1]
		alt_list = temp_alts.split(",")
		cov_list = cov.split(",")
		call_dict = {}
		for i in range(len(alt_list)):
			call_dict[alt_list[i]] = int(cov_list[i])
	
		total_cov = sum(call_dict.values())
		if total_cov<min_cov:
			base_calls[chrom][pos] = "NA"
			continue
		cutoff = 0
		if total_cov==min_cov:
			cutoff = min_cov
		elif total_cov>min_cov:
			cutoff = total_cov*read_frac
		allele = ""
		for call in call_dict:
			if call_dict[call]>=cutoff:
				allele = allele+call
		if allele=="":
			allele = "NA"
		if allele!=ref:
			base_calls[chrom][pos] = allele
		else:
			base_calls[chrom][pos] = "ref"
	return base_calls
	

def calls2variants(prefix,base_calls):
	
	dr_dict =load_bed(dr_bed_file)
	temp_bed = prefix+"_temp.bed"
	dr_calls = defaultdict(list)
	missing_calls = []
	with open(temp_bed,"w") as o:
		for chrom in sorted(base_calls):
			for pos in [str(y) for y in sorted([int(x) for x in base_calls[chrom]])]:
				if base_calls[chrom][pos]!="ref":
					obj = {"chr":chrom,"pos":pos,"alt":base_calls[chrom][pos],"drugs":dr_dict[chrom][pos]}
					if base_calls[chrom][pos]!="NA":
						o.write("%s\t%s\t%s\n" % (chrom,int(pos)-1,pos))
						for d in dr_dict[chrom][pos]:
							dr_calls[d].append(obj)
					elif base_calls[chrom][pos]=="NA":
						missing_calls.append(obj)
	return dr_calls,missing_calls


def annotate_calls(prefix,dr_calls):
	temp_bed = prefix+"_temp.bed"
	ann_dict = load_ann(temp_bed)
	final_results = []
	for drug in dr_calls:
		for var in dr_calls[drug]:
			temp_res = ""
			chrom = var["chr"]
			pos = var["pos"]
			gene = ann_dict[chrom][pos]["gene"]
			lt = ann_dict[chrom][pos]["rv"]
			temp_res += "%s\t%s\t%s\t%s\t%s\t%s\t" % (drug,prefix,chrom,pos,gene,lt)
			alt_nt = var["alt"]
			if len(alt_nt)>1:
				temp_res += "INDEL\t%s" % var["alt"]
				continue ## indel is found
			alt_aa = ann_dict[chrom][pos]["ann"][alt_nt]["aa"]	
			ref_aa = ann_dict[chrom][pos]["ref_aa"]
			codon_pos = ann_dict[chrom][pos]["codon_num"]	
			change = ref_aa+codon_pos+alt_aa
			if ref_aa==alt_aa:
				ref_nt = ann_dict[chrom][pos]["ref_nt"]
				gene_pos = ann_dict[chrom][pos]["gene_nt"]
				change=ref_nt+gene_pos+alt_nt
			temp_res += "SNP\t%s" % change
			final_results.append(temp_res)
	print "\n".join(final_results)
	open(prefix+".results.txt","w").write("\n".join(final_results))
def cleanup(prefix):
	cmd_clean = "rm %s.bam* %s.pileup* %s_temp.bed" % (prefix,prefix,prefix)
	subprocess.call(cmd_clean,shell=True)

################ Main functions ##################

def main_run_pipeline(args):
	check_files(args.read1,args.read2)
	doMapping(args.read1,args.read2,args.threads,args.prefix)		
	doPileup(args.prefix)
	base_calls = getCalls(prefix,args.min_cov,args.read_frac)
	dr_calls,missing_calls = calls2variants(args.prefix,base_calls)	
	annotate_calls(args.prefix,dr_calls)
	print missing_calls
	if args.clean:
		cleanup(args.prefix)
####################### Argument Parser ###########################

parser = argparse.ArgumentParser(description='Python wrapper to filter variants',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")


parser_sub = subparsers.add_parser('sv', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--read1','-1',help='First read file [required]')
parser_sub.add_argument('--read2','-2',help='Second read file')
parser_sub.add_argument('--prefix','-p',default="tbprofiler",help='Sample prefix')
parser_sub.add_argument('--min_cov','-c',type=int,default=10,help='Minimum coverage required')
parser_sub.add_argument('--read_frac','-f',type=float,default=0.7,help='Minimum fraction of total reads to call allele')
parser_sub.add_argument('--threads','-t',default="1",help='Threads')
parser_sub.add_argument('--clean',help="Remove temporary files",action="store_true")
parser_sub.set_defaults(func=main_run_pipeline)

#parser_sub = subparsers.add_parser('stats', help='Generate raw unfiltered matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#parser_sub.set_defaults(func=main_filter)


args = parser.parse_args()
args.func(args)















