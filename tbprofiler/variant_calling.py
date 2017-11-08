from __future__ import division
from files import *
from collections import defaultdict

min_dp = 10
min_frac = 0.7


def do_pileup(self,bed_file=None):
    self.params["temp"] = bed_file
    if bed_file:
        cmd = "%(samtools)s view -bL %(temp)s %(bamfile)s  > %(temp_bam)s && %(htsbox)s pileup -f %(reffile)s -Q 8 %(temp_bam)s > %(temp_pileup)s" % self.params
    else:
        cmd = "%(samtools)s view -b %(bamfile)s  > %(temp_bam)s && %(htsbox)s pileup -f %(reffile)s -Q 8 %(temp_bam)s > %(temp_pileup)s" % self.params
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
            #Chromosome      23      G       G,G-1C,G+3AAA   0/1:49,1,1
            arr = l.rstrip().split()
            if (arr[0],arr[1]) not in bed_pos: continue
            alleles = arr[3].split(",")
            depth = [int(x) for x in arr[4].split(":")[1].split(",")]
            max_allele_dp = max(depth)
            max_allele = alleles[depth.index(max_allele_dp)]
            max_allele_frac = max_allele_dp/sum(depth)
            if sum(depth)<min_dp:
                call = "N"
            if max_allele_frac<min_frac:
                call = "N"
            else:
                call = max_allele
            final_calls[arr[0]][arr[1]].append((call,1,max_allele_dp))
    return final_calls

def call_dr_variants(self):
    if self.params["platform"] == "Illumina":
        cmd = "set -euf pipefail; %(samtools)s view -bL %(dr_bed_file)s %(bamfile)s | %(lofreq)s indelqual -u 30 - |%(lofreq)s call -f %(reffile)s --call-indels --no-default-filter - | %(bcftools)s norm -f %(reffile)s -  | %(bcftools)s csq -f %(reffile)s -g %(gfffile)s - > %(dr_vcffile)s" %  self.params
        run_cmd(cmd,verbose=self.params["verbose"])
    if self.params["platform"] == "minION":
        pileup2vcf(self,bed_file=self.params["dr_bed_file"])
        cmd = "%(bcftools)s csq -g %(gfffile)s -f %(reffile)s %(temp_file)s > %(dr_vcffile)s" % self.params
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

    do_pileup(self,bed_file)

    OUT = open("%(temp_file)s" % self.params,"w")
    OUT.write(header)
    for l in open(self.params["temp_pileup"]):
        #Chromosome      23      G       G,G-1C,G+3AAA   0/1:49,1,1
        arr = l.rstrip().split()
        alleles = arr[3].split(",")
        depth = [int(x) for x in arr[4].split(":")[1].split(",")]

        max_allele_dp = max(depth)
        max_allele = alleles[depth.index(max_allele_dp)]
        max_allele_frac = max_allele_dp/sum(depth)
        DP4 = "0,%s,0,%s" % (sum(depth)-max_allele_dp,max_allele_dp)
        if sum(depth)<min_dp: continue
        if max_allele_frac<min_frac:
            call = "N"; continue
        else:
            call = max_allele
        if call==arr[2]: continue

        #Chromosome    7    .    G    A    8    .    .    GT:AD    0/1:914,8
        OUT.write("Chromosome\t%s\t.\t%s\t%s\t255\t.\tDP4=%s\tGT\t1\n" % (arr[1],arr[2],max_allele,DP4))
    OUT.close()
