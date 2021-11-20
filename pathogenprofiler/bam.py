from .utils import add_arguments_to_self, run_cmd, cmd_out, filecheck, index_bam
from .vcf import vcf, delly_bcf
from collections import defaultdict
import json
from uuid import uuid4
import os

class bam:
    """
    A class to perform operations on BAM files such as SNP calling
    """
    def __init__(self,bam_file,prefix,platform,threads=1):
        add_arguments_to_self(self, locals())
        filecheck(self.bam_file)
        index_bam(bam_file,threads=threads)

    def run_delly(self):
        if self.platform=="illumina":
            _,stderr = run_cmd("delly call -t DEL -g %(ref_file)s %(bam_file)s -o %(prefix)s.delly.bcf" % vars(self),terminate_on_error=False)
            if "not enough data to estimate library parameters" in stderr:
                return None
            else:
                return delly_bcf("%(prefix)s.delly.bcf" % vars(self))
        else:
            _,stderr = run_cmd("delly lr -t DEL -g %(ref_file)s %(bam_file)s -o %(prefix)s.delly.bcf" % vars(self),terminate_on_error=False)
            if "not enough data to estimate library parameters" in stderr:
                return None
            else:
                return delly_bcf("%(prefix)s.delly.bcf" % vars(self))
                
    def call_variants(self,ref_file,caller,bed_file=None,threads=1,calling_params=None,remove_missing=False, samclip=False,min_dp=10):
        add_arguments_to_self(self, locals())
        filecheck(ref_file)
        self.caller = caller.lower()
        self.missing_cmd = "" if remove_missing else "-S ."

        # Set up final vcf file name
        # Make the windows for parallel calling based on chunking the whole
        # genome or by providing a bed file
        if bed_file:
            self.vcf_file = "%s.targets.vcf.gz" % (self.prefix) if bed_file else "%s.vcf.gz" % (self.prefix)
            self.windows_cmd = "cat %(bed_file)s | awk '{print $1\":\"$2\"-\"$3\" \"$1\"_\"$2\"_\"$3}'" % vars(self)
        else:
            self.vcf_file = "%s.vcf.gz" % (self.prefix) if bed_file else "%s.vcf.gz" % (self.prefix)
            self.windows_cmd = "bedtools makewindows -g %(ref_file)s.fai -n %(threads)s | awk '{print $1\":\"$2+1\"-\"$3\" \"$1\"_\"$2+1\"_\"$3}'" % vars(self)

        self.samclip_cmd = "| samclip --ref %(ref_file)s" % vars(self) if samclip else ""
        # Run through different options. Start with nanopore because it should
        # only be run with bcftools and will not even take caller into account
        # if it is set the the wrong option
        if self.platform == "nanopore" and self.caller=="bcftools":
            self.calling_params = calling_params if calling_params else "-Bq8"
            self.calling_cmd = "bcftools mpileup -f %(ref_file)s %(calling_params)s -a DP,AD -r {1} %(bam_file)s | bcftools call -mv | bcftools view -c 1 | bcftools norm -f %(ref_file)s | bcftools filter -e 'FMT/DP<%(min_dp)s' %(missing_cmd)s | bcftools filter -e 'IMF < 0.7' -S 0 -Oz -o %(prefix)s.{2}.vcf.gz" % vars(self)
        elif self.platform == "nanopore" and self.caller=="freebayes":
            self.calling_params = calling_params if calling_params else "-F 0.7"
            self.calling_cmd = "freebayes -f %(ref_file)s -r {1} --haplotype-length -1 %(calling_params)s %(bam_file)s | bcftools view -c 1 | bcftools norm -f %(ref_file)s | bcftools filter -t {1} -e 'FMT/DP<%(min_dp)s' %(missing_cmd)s -Oz -o %(prefix)s.{2}.vcf.gz" % vars(self)
        elif self.platform=="illumina" and self.caller == "bcftools":
            self.calling_params = calling_params if calling_params else "-ABq8 -Q0"
            self.calling_cmd = "samtools view -T %(ref_file)s -h %(bam_file)s {1} %(samclip_cmd)s | samtools view -b > %(prefix)s.{2}.tmp.bam && samtools index %(prefix)s.{2}.tmp.bam && bcftools mpileup -f %(ref_file)s %(calling_params)s -a DP,AD -r {1} %(prefix)s.{2}.tmp.bam | bcftools call -mv | bcftools norm -f %(ref_file)s | bcftools filter -t {1} -e 'FMT/DP<%(min_dp)s' %(missing_cmd)s -Oz -o %(prefix)s.{2}.vcf.gz" % vars(self)
        elif self.platform=="illumina" and self.caller == "gatk":
            self.calling_params = calling_params if calling_params else ""
            self.calling_cmd = "samtools view -T %(ref_file)s -h %(bam_file)s {1} %(samclip_cmd)s | samtools view -b > %(prefix)s.{2}.tmp.bam && samtools index %(prefix)s.{2}.tmp.bam && gatk HaplotypeCaller -R %(ref_file)s -I %(prefix)s.{2}.tmp.bam -O /dev/stdout -L {1} %(calling_params)s -OVI false | bcftools norm -f %(ref_file)s | bcftools filter -t {1} -e 'FMT/DP<%(min_dp)s' %(missing_cmd)s -Oz -o %(prefix)s.{2}.vcf.gz" % vars(self)
        elif self.platform=="illumina" and self.caller == "freebayes":
            self.calling_params = calling_params if calling_params else ""
            self.calling_cmd = "samtools view -T %(ref_file)s -h %(bam_file)s {1} %(samclip_cmd)s | samtools view -b > %(prefix)s.{2}.tmp.bam && samtools index %(prefix)s.{2}.tmp.bam && freebayes -f %(ref_file)s -r {1} --haplotype-length -1 %(calling_params)s %(prefix)s.{2}.tmp.bam | bcftools view -c 1 | bcftools norm -f %(ref_file)s | bcftools filter -t {1} -e 'FMT/DP<%(min_dp)s' %(missing_cmd)s -Oz -o %(prefix)s.{2}.vcf.gz" % vars(self)


        run_cmd('%(windows_cmd)s | parallel -j %(threads)s --col-sep " " "%(calling_cmd)s"' % vars(self))
        run_cmd('%(windows_cmd)s | parallel -j %(threads)s --col-sep " " "bcftools index  %(prefix)s.{2}.vcf.gz"' % vars(self) )
        run_cmd("bcftools concat -aD `%(windows_cmd)s | awk '{print \"%(prefix)s.\"$2\".vcf.gz\"}'` | bcftools view -c1 -a -Oz -o %(vcf_file)s" % vars(self))
        run_cmd("rm `%(windows_cmd)s | awk '{print \"%(prefix)s.\"$2\".vcf.gz*\"}'`" % vars(self))
        if self.platform!="nanopore":
            run_cmd("rm `%(windows_cmd)s | awk '{print \"%(prefix)s.\"$2\".tmp.bam*\"}'`" % vars(self))

        return vcf(self.vcf_file)

    def flagstat(self):
        tmpfile = str(uuid4())
        run_cmd(f"samtools flagstat -O json {self.bam_file} > {tmpfile}")
        flagstat = json.load(open(tmpfile))
        self.num_reads_mapped = flagstat["QC-passed reads"]["mapped"]
        self.pct_reads_mapped = round(flagstat["QC-passed reads"]["mapped"]/flagstat["QC-passed reads"]["total"]*100,2)
        os.remove(tmpfile)
        return self.num_reads_mapped,self.pct_reads_mapped

    def get_median_coverage(self):
        lines = []
        for l in cmd_out("bedtools genomecov -ibam %s" % (self.bam_file)):
            arr = l.split()
            if arr[0]=="genome":
                lines.append(arr)
        midpoint =  int(lines[0][3])/2
        x = 0
        for row in lines:
            x = x + int(row[2])
            if x>midpoint:
                break
        self.median_coverage = int(row[1])
        return int(row[1])

    def bed_zero_cov_regions(self,bed_file):
        add_arguments_to_self(self, locals())
        cmd = "bedtools coverage -sorted -b %(bam_file)s -a %(bed_file)s" % vars(self)
        results = []
        for l in cmd_out(cmd):
            results.append(l.split())
        return results

    def get_bed_gt(self,bed_file,ref_file,caller,platform):
        add_arguments_to_self(self, locals())
        results = defaultdict(lambda : defaultdict(dict))
        run_cmd("samtools view -Mb -L %(bed_file)s %(bam_file)s > %(prefix)s.tmp.bam" % vars(self))
        run_cmd("samtools index %(prefix)s.tmp.bam" % vars(self))
        if platform=="nanopore":
            caller="bcftools"
        if caller == "gatk":
            cmd = "gatk HaplotypeCaller -I %(prefix)s.tmp.bam -R %(ref_file)s -L %(bed_file)s -OVI false -O /dev/stdout | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD]\\n'" % vars(self)
        elif caller == "freebayes":
            cmd = "freebayes -f %(ref_file)s -t %(bed_file)s %(prefix)s.tmp.bam --haplotype-length -1 | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD]\\n'" % vars(self)
        else:
            cmd = "bcftools mpileup -f %(ref_file)s -R %(bed_file)s %(prefix)s.tmp.bam -BI -a AD | bcftools call -mv | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD]\\n'" % vars(self)

        for l in cmd_out(cmd):
            # Chromosome    4348079    0/0    51
            chrom, pos, ref, alt, gt, ad = l.rstrip().split()
            pos = int(pos)
            d = {}
            alts = alt.split(",")
            ad = [int(x) for x in ad.split(",")]
            if gt == "0/0":
                d[ref] = ad[0]
            elif gt == "./.":
                d[ref] = 0
            else:
                genotypes = list([ref]+alts)
                if self.platform == "nanopore":
                    idx = ad.index(max(ad))
                    d[genotypes[idx]] = ad[idx]
                else:
                    for i, a in enumerate(genotypes):
                        d[a] = ad[i]
            results[chrom][pos] = d

        ref_nt = {}
        for l in cmd_out("bedtools getfasta -fi %s -bed %s" % (ref_file,bed_file)):
            if l[0]==">":
                tmp = l.strip().replace(">","").split(":")
                tmp_chrom = tmp[0]
                tmp_pos = int(tmp[1].split("-")[1])
            else:
                ref_nt[(tmp_chrom,tmp_pos)] = l.strip()

        for l in cmd_out(f"samtools view -Mb -L {bed_file} {self.bam_file} | bedtools coverage -a {bed_file} -b - -d -sorted"):
            row = l.strip().split()
            chrom = row[0]
            pos = int(row[2])
            cov = int(row[-1])
            if chrom not in results or pos not in results[chrom]:
                results[chrom][pos] = {ref_nt[(chrom,pos)]:cov}

        return results

    def get_region_coverage(self,bed_file,per_base=False,group_column=4,fraction_threshold=0,region_column=3):
        add_arguments_to_self(self, locals())
        self.region_cov = defaultdict(list)
        self.region_fraction = []
        self.genome_coverage = []

        for l in cmd_out(f"samtools view -Mb -L {bed_file} {self.bam_file} | bedtools coverage -a {bed_file} -b - -d -sorted"):

            row = l.split()
            region = row[region_column]
            depth = int(row[-1])
            genomic_position = int(row[1]) + int(row[-2]) -1
            self.genome_coverage.append((genomic_position, depth))
            self.region_cov[region].append(depth)

        for region in self.region_cov:
            total_num = len(self.region_cov[region])
            num_thresh = len([d for d in self.region_cov[region] if d<=fraction_threshold])
            self.region_fraction.append({"gene_id":region, "fraction":num_thresh/total_num, "cutoff": fraction_threshold})
        if per_base:
            return self.region_cov
        else:
            return self.region_fraction

    def get_missing_genomic_positions(self,bed_file=None,cutoff=10):
        if not hasattr(self,"genome_coverage"):
            self.get_region_coverage(bed_file)
        return [x[0] for x in self.genome_coverage if x[1]<cutoff]
