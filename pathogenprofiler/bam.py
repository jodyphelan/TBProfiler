from .utils import add_arguments_to_self, run_cmd, cmd_out, filecheck, index_bam, nofile, rm_files, log, load_bed, median
from .fasta import fasta
from .vcf import vcf, delly_bcf
from tqdm import tqdm
from collections import defaultdict
import sys
class bam:
    """
    A class to perform operations on BAM files such as SNP calling
    """
    def __init__(self,bam_file,prefix,platform,threads=1):
        add_arguments_to_self(self, locals())
        filecheck(self.bam_file)
        index_bam(bam_file,threads=threads)

    def run_delly(self):
        stdout,stderr = run_cmd("delly call -t DEL -g %(ref_file)s %(bam_file)s -o %(prefix)s.delly.bcf" % vars(self),terminate_on_error=False)
        if "not enough data to estimate library parameters" in stderr:
            return None
        else:
            return delly_bcf("%(prefix)s.delly.bcf" % vars(self))

    def call_variants(self,ref_file,caller,bed_file=None,threads=1,calling_params=None,remove_missing=False,min_depth=10):
        add_arguments_to_self(self, locals())
        filecheck(ref_file)
        self.caller = caller.lower()
        self.missing_cmd = "" if remove_missing else "-S ."
        # Set up final vcf file name
        self.vcf_file = "%s.targets.vcf.gz" % (self.prefix) if bed_file else "%s.vcf.gz" % (self.prefix)

        # Make the windows for parallel calling based on chunking the whole
        # genome or by providing a bed file
        if bed_file:
            self.windows_cmd = "cat %(bed_file)s | awk '{print $1\":\"$2\"-\"$3\" \"$1\"_\"$2\"_\"$3}'" % vars(self)
        else:
            self.windows_cmd = "bedtools makewindows -g %(ref_file)s.fai -n %(threads)s | awk '{print $1\":\"$2+1\"-\"$3\" \"$1\"_\"$2+1\"_\"$3}'" % vars(self)

        # Run through different options. Start with nanopore because it should
        # only be run with bcftools and will not even take caller into account
        # if it is set the the wrong option
        if self.platform == "nanopore":
            self.calling_params = calling_params if calling_params else "-Bq8"
            self.calling_cmd = "bcftools mpileup -f %(ref_file)s %(calling_params)s -a DP,AD -r {1} %(bam_file)s | bcftools call -mv | bcftools filter -e 'FMT/DP<%(min_depth)s' %(missing_cmd)s | bcftools filter -e 'IMF < 0.7' -S 0 -Oz -o %(prefix)s.{2}.vcf.gz" % vars(self)
        elif self.caller == "bcftools":
            self.calling_params = calling_params if calling_params else "-ABq8 -Q0"
            self.calling_cmd = "bcftools mpileup -f %(ref_file)s %(calling_params)s -a DP,AD -r {1} %(bam_file)s | bcftools call -mv | bcftools norm -f %(ref_file)s | bcftools filter -e 'FMT/DP<%(min_depth)s' %(missing_cmd)s -Oz -o %(prefix)s.{2}.vcf.gz" % vars(self)
        elif self.caller == "gatk":
            self.calling_params = calling_params if calling_params else ""
            self.calling_cmd = "gatk HaplotypeCaller -R %(ref_file)s -I %(bam_file)s -O %(prefix)s.{2}.vcf.gz -L {1} %(calling_params)s" % vars(self)
        elif self.caller == "freebayes":
            self.calling_params = calling_params if calling_params else ""
            self.calling_cmd = "freebayes -f %(ref_file)s -r {1} %(bam_file)s --haplotype-length 1 %(calling_params)s | bcftools view -c 1 | bcftools norm -f %(ref_file)s | bcftools filter -e 'FMT/DP<%(min_depth)s' %(missing_cmd)s -Oz -o %(prefix)s.{2}.vcf.gz" % vars(self)

        run_cmd('%(windows_cmd)s | parallel -j %(threads)s --col-sep " " "%(calling_cmd)s"' % vars(self))
        run_cmd('%(windows_cmd)s | parallel -j %(threads)s --col-sep " " "bcftools index  %(prefix)s.{2}.vcf.gz"' % vars(self) )
        run_cmd("bcftools concat -aD -Oz -o %(vcf_file)s `%(windows_cmd)s | awk '{print \"%(prefix)s.\"$2\".vcf.gz\"}'`" % vars(self))
        run_cmd("rm `%(windows_cmd)s | awk '{print \"%(prefix)s.\"$2\".vcf.gz*\"}'`" % vars(self))

        return vcf(self.vcf_file)

    def flagstat(self):
        lines = []
        for l in cmd_out("samtools flagstat %s" % (self.bam_file)):
            arr = l.split()
            lines.append(arr)
        self.num_reads_mapped = int(lines[4][0])
        self.pct_reads_mapped = 0.0 if self.num_reads_mapped == 0 else float(lines[4][4][1:-1])
        return self.num_reads_mapped,self.pct_reads_mapped

    def bed_zero_cov_regions(self,bed_file):
        add_arguments_to_self(self, locals())
        cmd = "bedtools coverage -b %(bam_file)s -a %(bed_file)s" % vars(self)
        results = []
        for l in cmd_out(cmd):
            results.append(l.split())
        return results

    def get_bed_gt(self,bed_file,ref_file,caller):
        add_arguments_to_self(self, locals())
        results = defaultdict(lambda : defaultdict(dict))

        if caller == "gatk":
            cmd = "gatk HaplotypeCaller -I %(bam_file)s -R %(ref_file)s -L %(bed_file)s -ERC BP_RESOLUTION -OVI false -O /dev/stdout | bcftools view -a | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD]\\n'" % vars(self)
        else:
            cmd = "bcftools mpileup -f %(ref_file)s -R %(bed_file)s %(bam_file)s -BI -a AD | bcftools call -m | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD]\\n'" % vars(self)
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

        for l in cmd_out("bedtools coverage -a %s -b %s -d -sorted" % (bed_file,self.bam_file)):
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
        for l in cmd_out("bedtools coverage -a %(bed_file)s -b %(bam_file)s -d -sorted" % vars(self)):
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

    def get_missing_genomic_positions(self,cutoff=10):
        if not hasattr(self,"genome_coverage"):
            self.get_region_coverage()
        return [x[0] for x in self.genome_coverage if x[1]<cutoff]
