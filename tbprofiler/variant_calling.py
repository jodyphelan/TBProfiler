from __future__ import division
from files import *
from collections import defaultdict
import re
min_dp = 10
min_frac = 0.7

indelre = re.compile("(\w)[\+\-](\d+)(\w+)")

def recode_indels(indels):
	#["C+5CGGGG","G-1C"]
	sorted_indels = sorted([x for x in indels],key=lambda y:len(y))
	largest_del = sorted([x for x in indels if "-" in x],key=lambda y:len(y))

	if len(largest_del)==0:
		leftseq = indelre.search(sorted_indels[-1]).group(1)
		rightseq = ""
	else:
		leftseq = indelre.search(largest_del[-1]).group(1)
		rightseq = indelre.search(largest_del[-1]).group(3)
	refseq = leftseq+rightseq
	recoded_indels = []
	for i in indels:
		if "-" in i:
			indel_len = int(indelre.search(i).group(2))
			iseq = leftseq+rightseq[:-(len(rightseq)-indel_len)]
		elif "+" in i:
			iseq = leftseq+indelre.search(i).group(3)+rightseq
		else:
			iseq = i+rightseq
		recoded_indels.append(iseq)
	return (refseq,recoded_indels)

def do_pileup(self,bed_file=None):
	self.params["temp"] = bed_file
	if bed_file:
		cmd = "%(samtools)s view -@ %(threads)s -bL %(temp)s %(bamfile)s  > %(temp_bam)s && %(htsbox)s pileup -f %(reffile)s -Q 8 %(temp_bam)s > %(temp_pileup)s" % self.params
	else:
		cmd = "%(htsbox)s pileup -f %(reffile)s -Q 8 %(bamfile)s > %(temp_pileup)s" % self.params
	run_cmd(cmd,verbose=self.params["verbose"])

def htsbox_calls(self,bed_file):
	bed_pos = set()
	for l in open(bed_file):
		arr = l.rstrip().split()
		for i in range(int(arr[1]),int(arr[2])+1):
			bed_pos.add((arr[0],str(i)))
	do_pileup(self,bed_file=bed_file)
	final_calls = defaultdict(lambda :defaultdict(list))
	if self.params["platform"] == "Illumina":
		for l in open(self.params["temp_pileup"]):
			arr = l.rstrip().split()
			if (arr[0],arr[1]) not in bed_pos: continue
			calls = arr[3].split(",")
			cov = [int(x) for x in arr[4].split(":")[1].split(",")]
			tot = sum(cov)
			for i in range(len(calls)):
				final_calls[arr[0]][arr[1]].append((calls[i],cov[i]/tot,cov[i]))
	elif self.params["platform"] == "minION":
		for l in open(self.params["temp_pileup"]):
			#Chromosome	  23	  G	   G,G-1C,G+3AAA   0/1:49,1,1
			arr = l.rstrip().split()
			if (arr[0],arr[1]) not in bed_pos: continue
			alleles = arr[3].split(",")
			depth = [int(x) for x in arr[4].split(":")[1].split(",")]
			max_allele_dp = max(depth)
			max_allele = alleles[depth.index(max_allele_dp)]
			max_allele_frac = max_allele_dp/sum(depth)
			if len(max_allele)>1:
				max_allele = recode_indels([arr[1],max_allele])[1][0]
			if sum(depth)<min_dp:
				call = "N"
			if max_allele_frac<min_frac:
				call = "N"
			else:
				call = max_allele
			final_calls[arr[0]][arr[1]].append((call,1,max_allele_dp))
	return final_calls

def call_variants(self,bed_file=False,vcf_file=False,caller="bcftools",gvcf=False):
	if bed_file:
		self.params["temp_bed"] = bed_file
		self.params["temp_vcf_file"] = vcf_file
	if gvcf: caller="bcftools"

	if self.params["platform"] == "Illumina":
		if caller=="lofreq":
			if bed_file:
				cmd = "set -euf pipefail; %(samtools)s view -bL %(temp_bed)s %(bamfile)s | %(lofreq)s indelqual -u 30 - |%(lofreq)s call -f %(reffile)s --call-indels --no-default-filter - | bcftools view -i 'AF>%(af)s' | %(bcftools)s norm -f %(reffile)s - | %(bcftools)s csq -f %(reffile)s -g %(gfffile)s - > %(temp_vcf_file)s" %  self.params
			else:
				cmd = "set -euf pipefail; %(samtools)s view -b %(bamfile)s | %(lofreq)s indelqual -u 30 - |%(lofreq)s call -f %(reffile)s --call-indels --no-default-filter - | bcftools view -i 'AF>%(af)s' | %(bcftools)s norm -f %(reffile)s -  | %(bcftools)s csq -f %(reffile)s -g %(gfffile)s - > %(vcffile)s" %  self.params
		elif caller=="bcftools":
			if bed_file:
				cmd = "set -euf pipefail; %(samtools)s view -bL %(temp_bed)s %(bamfile)s | %(samtools)s mpileup -ugf %(reffile)s -t DP,AD - | %(bcftools)s call --threads %(threads)s -mv  | bcftools view -i 'AF>%(af)s' | %(bcftools)s norm -f %(reffile)s -  | %(bcftools)s csq --phase a -f %(reffile)s -g %(gfffile)s - > %(temp_vcf_file)s" %  self.params
			else:
				if gvcf:
					cmd = "set -euf pipefail; %(samtools)s view -b %(bamfile)s | %(samtools)s mpileup -ugf %(reffile)s -t DP,AD - | %(bcftools)s call --threads %(threads)s -mg 10 | bcftools view -i 'AF>%(af)s' | %(bcftools)s norm -f %(reffile)s -O z -o %(gvcffile)s" %  self.params
				else:
					cmd = "set -euf pipefail; %(samtools)s view -b %(bamfile)s | %(samtools)s mpileup -ugf %(reffile)s -t DP,AD - | %(bcftools)s call --threads %(threads)s -mv | bcftools view -i 'AF>%(af)s' | %(bcftools)s norm -f %(reffile)s -  | %(bcftools)s csq --phase a -f %(reffile)s -g %(gfffile)s - > %(vcffile)s" %  self.params
		else:
			print("\nOnly lofreq or bcftools available to call variants...exiting\n"); quit()

		run_cmd(cmd,verbose=self.params["verbose"])
	elif self.params["platform"] == "minION":
		if bed_file:
			pileup2vcf(self,bed_file=self.params["temp_bed"])
			cmd = "%(bcftools)s csq -g %(gfffile)s -f %(reffile)s %(temp_file)s > %(temp_vcf_file)s" % self.params
		else:
			pileup2vcf(self)
			cmd = "%(bcftools)s csq -g %(gfffile)s -f %(reffile)s %(temp_file)s > %(vcffile)s" % self.params
		run_cmd(cmd,verbose=self.params["verbose"])


def pileup2vcf(self,bed_file=False):
	header = """##fileformat=VCFv4.1
##source=htsbox-pileup-r340
##reference=/home/jody/refgenome/MTB-h37rv_asm19595v2-eg18.fa
##contig=<ID=Chromosome,length=4411532>
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%(prefix)s
""" % self.params

	if bed_file:
		do_pileup(self,bed_file)
		bed_pos = set()
		for l in open(bed_file):
			arr = l.rstrip().split()
			for i in range(int(arr[1]),int(arr[2])):
				bed_pos.add(str(i))
	else:
		do_pileup(self)
	OUT = open("%(temp_file)s" % self.params,"w")
	OUT.write(header)
	for l in open(self.params["temp_pileup"]):
		#Chromosome	  23	  G	   G,G-1C,G+3AAA   0/1:49,1,1
		arr = l.rstrip().split()
		if bed_file:
			if arr[1] not in bed_pos:
				continue
		alleles = arr[3].split(",")
		depth = [int(x) for x in arr[4].split(":")[1].split(",")]
		ref = arr[2]
		max_allele_dp = max(depth)
		max_allele = alleles[depth.index(max_allele_dp)]
		max_allele_frac = max_allele_dp/sum(depth)
		adjusted_allele_frac = max_allele_dp/(max_allele_dp+sorted(depth)[-2]) if len(depth)>1 else max_allele_frac
		if len(max_allele)>1:
			indel = recode_indels([max_allele])

			max_allele = indel[1][0]
			ref = indel[0]

		DP4 = "0,%s,0,%s" % (sum(depth)-max_allele_dp,max_allele_dp)
		if sum(depth)<min_dp: continue

#		if max_allele_frac<min_frac:
		if adjusted_allele_frac<min_frac:
			call = "N"; continue
		else:
			call = max_allele

		if call==ref: continue

		#Chromosome	7	.	G	A	8	.	.	GT:AD	0/1:914,8
		OUT.write("Chromosome\t%s\t.\t%s\t%s\t255\t.\tDP4=%s\tGT\t1\n" % (arr[1],ref,max_allele,DP4))
	OUT.close()
