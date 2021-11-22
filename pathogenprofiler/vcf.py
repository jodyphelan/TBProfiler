from .utils import run_cmd, cmd_out,add_arguments_to_self,rm_files, index_bcf,tabix, log, load_bed, debug, warninglog
from .fasta import fasta
from collections import defaultdict
import re
from uuid import uuid4
import sys

re_seq = re.compile("([0-9\-]*)([A-Z\*]+)")
re_I = re.compile("([A-Z\*]+)")
number_re = re.compile("[0-9\-]+")

def parse_mutation(x):
    tmp = x.split(">")
    aa_changed = True if len(tmp)>1 else False
    re_obj = re_seq.search(tmp[0])
    change_num = re_obj.group(1)
    ref_aa = re_obj.group(2)
    alt_aa = re_seq.search(tmp[1]).group(2) if aa_changed else None
    return change_num,ref_aa,alt_aa

class vcf:
    def __init__(self,filename,prefix=None,threads=1):
        self.samples = []
        add_arguments_to_self(self,locals())
        if prefix==None:
            if filename[-4:] == ".bcf":
                self.prefix = filename[:-4]
            elif filename[-7:] == ".vcf.gz":
                self.prefix = filename[:-7]
            elif filename[-8:] == ".gvcf.gz":
                self.prefix = filename[:-8]
            elif filename[-4:] == ".vcf":
                self.prefix = filename[:-4]
            else:
                self.prefix = filename
        else:
            self.prefix = prefix
        if filename[-4:] == ".bcf":
            index_bcf(filename,threads)
        else:
            tabix(filename,threads)
        for l in cmd_out("bcftools query -l %(filename)s" % vars(self)):
            self.samples.append(l.rstrip())
    def view_regions(self,bed_file):
        self.bed_file = bed_file
        self.newfile = "%(prefix)s.targets.vcf.gz" % vars(self)
        run_cmd("bcftools view -R %(bed_file)s %(filename)s -Oz -o %(newfile)s" % vars(self))
        return vcf(self.newfile)


    def run_snpeff(self,db,ref_file,gff_file,rename_chroms = None, split_indels=True):
        add_arguments_to_self(self,locals())
        self.vcf_csq_file = self.prefix+".csq.vcf.gz"
        self.rename_cmd = f"rename_vcf_chrom.py --source {' '.join(rename_chroms['source'])} --target {' '.join(rename_chroms['target'])} |" if rename_chroms else ""
        self.re_rename_cmd = f"| rename_vcf_chrom.py --source {' '.join(rename_chroms['target'])} --target {' '.join(rename_chroms['source'])}" if rename_chroms else ""
        if split_indels:
            self.tmp_file1 = "%s.vcf" % uuid4()
            self.tmp_file2 = "%s.vcf" % uuid4()

            run_cmd("bcftools view -v snps %(filename)s | combine_vcf_variants.py --ref %(ref_file)s --gff %(gff_file)s | %(rename_cmd)s snpEff -Djava.net.useSystemProxies=true ann -noStats %(db)s - %(re_rename_cmd)s > %(tmp_file1)s" % vars(self))
            run_cmd("bcftools view -v indels %(filename)s | %(rename_cmd)s snpEff -Djava.net.useSystemProxies=true ann -noStats %(db)s - %(re_rename_cmd)s > %(tmp_file2)s" % vars(self))
            run_cmd("bcftools concat %(tmp_file1)s %(tmp_file2)s | bcftools sort -Oz -o %(vcf_csq_file)s" % vars(self))
            rm_files([self.tmp_file1, self.tmp_file2])


        else :
            run_cmd("bcftools view %(filename)s | %(rename_cmd)s snpEff -Djava.net.useSystemProxies=true ann %(db)s - %(re_rename_cmd)s | bcftools view -Oz -o %(vcf_csq_file)s" % vars(self))
        return vcf(self.vcf_csq_file,self.prefix)

        # run_cmd(f"snpEff ann {db} {self.filename} | bcftools view -Oz -o {self.prefix}.ann.vcf.gz")
        # return vcf(f"{self.prefix}.ann.vcf.gz")


    def load_ann(self,max_promoter_length=200, bed_file=None,intergenic=False,intragenic=False,upstream=False,downstream=False,noncoding=False,intronic=False,synonymous=False,splice=False,ablation=False):
        filter_out = []
        if intergenic==False:
            filter_out.append("intergenic_region")
        if intragenic==False:
            filter_out.append("intragenic_variant")
        if noncoding==False:
            filter_out.append("non_coding_transcript_variant")
            filter_out.append("non_coding_transcript_exon_variant")
        if downstream==False:
            filter_out.append("downstream_gene_variant")
        if upstream==False:
            filter_out.append("upstream_gene_variant")
        if intronic==False:
            filter_out.append("intron_variant")
        if synonymous==False:
            filter_out.append("synonymous_variant")
        if splice==False:
            filter_out.append("splice_region_variant&intron_variant")
        if ablation==False:
            filter_out.append("transcript_ablation")
        
        if bed_file:
            genes_to_keep = []
            for l in open(bed_file):
                row = l.strip().split()
                genes_to_keep.append((row[3],row[4]))

        variants = []
        for l in cmd_out(f"bcftools query -u -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%ANN\\t[%AD]\\n' {self.filename}"):
            chrom,pos,ref,alt_str,ann_str,ad_str = l.strip().split()

            alleles = [ref] + alt_str.split(",")
            if alt_str=="<DEL>":
                af_dict = {"<DEL>":1.0}
            else:
                ad = [int(x) for x in ad_str.split(",")]
                af_dict = {alleles[i]:ad[i]/sum(ad) for i in range(len(alleles))}
            ann_list = [x.split("|") for x in ann_str.split(",")]
            for alt in alleles[1:]:
                tmp_var = {
                    "chrom": chrom,
                    "genome_pos": int(pos),
                    "ref": ref,
                    "alt":alt,
                    "freq":af_dict[alt],
                    "consequences":[]
                }
                # if pos=="1473246":
                #     import pdb; pdb.set_trace()
                for ann in ann_list:
                    if ann[0]!=alt:
                        continue
                    if ann[1] in filter_out:
                        continue
                    if bed_file:
                        if ann[3] in [x[1] for x in genes_to_keep] or ann[4] in [x[0] for x in genes_to_keep]:
                            pass
                        else:
                            continue
                    if ann[1]=="upstream_gene_variant":
                        r = re.search("[cn].-([0-9]+)",ann[9])
                        if int(r.group(1))>max_promoter_length:
                            continue

                    tmp = {
                        "gene_name":ann[3],
                        "gene_id":ann[4],
                        "feature_id":ann[6],
                        "type":ann[1],
                        "nucleotide_change":ann[9],
                        "protein_change":ann[10],
                    }
                    tmp_var["consequences"].append(tmp)
                variants.append(tmp_var)
                    
        return variants

    
    def add_annotations(self,ref_file,bam_file):
        add_arguments_to_self(self,locals())
        self.new_file = self.prefix + ".ann.vcf.gz"

        run_cmd("gatk VariantAnnotator -R %(ref_file)s -I %(bam_file)s -V %(filename)s -O %(new_file)s -A MappingQualityRankSumTest -A ReadPosRankSumTest -A QualByDepth -A BaseQualityRankSumTest -A TandemRepeat -A StrandOddsRatio -OVI false" % vars(self))
        return vcf(self.new_file,self.prefix)


    def load_variants(self):
        variants = defaultdict(lambda:defaultdict(lambda:defaultdict(dict)))
        raw_variants = defaultdict(lambda:defaultdict(lambda:defaultdict(dict)))
        cmd = "bcftools query -u -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%TGT:%%AD]\\n' %s  | sed 's/\.\/\./N\/N/g'" % self.filename
        for l in cmd_out(cmd):
            row = l.split()
            alts = row[3].split(",")
            alleles = [row[2]]+alts
            for i in range(len(self.samples)):
                calls,ad = row[i+4].split(":")
                if calls=="N/N":
                    raw_variants[row[0]][row[1]][self.samples[i]]["N"] = 1.0
                    continue
                elif calls=="%s/%s" % (row[2],row[2]) and ad==".":
                    raw_variants[row[0]][row[1]][self.samples[i]][row[2]] = 1.0
                    continue
                ad = [int(x) if x!="." else 0 for x in ad.split(",")] if ad!="." else [0,100]
                sum_ad = sum(ad)
                for j in range(1,len(alleles)):
                    if ad[j]==0: continue
                    raw_variants[row[0]][row[1]][self.samples[i]][alleles[j]] = ad[j]/sum_ad
        for tchrom in raw_variants:
            for tpos in raw_variants[tchrom]:
                variants[tchrom][int(tpos)] = raw_variants[tchrom][tpos]
        return variants

    def get_positions(self):
        results = []
        for l in cmd_out("bcftools query -f '%%CHROM\\t%%POS\\n' %s" % self.filename):
            row = l.split()
            results.append((row[0],int(row[1])))
        return results

    def get_bed_gt(self,bed_file,ref_file):
        add_arguments_to_self(self,locals())
        cmd = "bcftools convert --gvcf2vcf -f %(ref_file)s %(filename)s  | bcftools view -T %(bed_file)s  | bcftools query -u -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD]\\n'" % vars(self)
        results = defaultdict(lambda : defaultdict(dict))
        ref_seq = fasta(ref_file).fa_dict
        for l in cmd_out(cmd):
            #Chromosome    4348079    0/0    51
            chrom,pos,ref,alt,gt,ad = l.rstrip().split()
            pos =int(pos)
            d = {}
            alts = alt.split(",")
            ad = [int(x) for x in ad.split(",")] if ad!="." else [0,100]
            if gt=="0/0":
                d[ref] = ad[0]
            elif gt=="./.":
                d[ref] = 0
            else:
                for i,a in enumerate([ref]+alts):
                    d[a] = ad[i]
            results[chrom][pos] = d
        bed = load_bed(bed_file,[1,3,5],1,3)
        for chrom in bed:
            for pos in bed[chrom]:
                if int(pos) not in results[chrom]:
                    results[chrom][int(pos)] = {ref_seq[chrom][pos-1]:50}
        return results

    def get_gatk_annotations(self):
        possible_annotations = ["ReadPosRankSum","QD","MQRankSum","ClippingRankSum","BaseQRankSum","SOR","TandemRepeat"]
        lines = []
        for l in cmd_out("bcftools view -h %s | grep INFO" % self.filename):
            lines.append(l)
        lines = "".join(lines)
        found_annotations = []
        for x in possible_annotations:
            if x in lines:
                found_annotations.append(x)
        return found_annotations

class delly_bcf(vcf):
    def __init__(self,filename):
         vcf.__init__(self,filename)
    def get_robust_calls(self,prefix,bed_file = None):
        self.tmpfile = f"{prefix}.tmp.delly.vcf.gz"
        self.outfile = f"{prefix}.delly.vcf.gz"
        run_cmd("bcftools view -c 2  %(filename)s | bcftools view -e '(INFO/END-POS)>=100000' -Oz -o %(tmpfile)s" % vars(self))
        if bed_file:
            run_cmd(f"bcftools index {self.tmpfile}")
            run_cmd(f"bcftools view -R {bed_file} {self.tmpfile} -Oz -o {self.outfile}")
        else:
            self.outfile = self.tmpfile
        return delly_bcf(self.outfile)

    # def overlap_bed(self,bed_file):
    #     results = []
    #     bed = load_bed(bed_file,[1,2,3,4,5],4)
    #     calls = self.get_robust_calls()
    #     for call in calls:
    #         set_call_pos = set(range(int(call[1]),int(call[2])))
    #         for region in bed:
    #             if bed[region][0]!=call[0]: continue
    #             set_region_pos = set(range(int(bed[region][1]),int(bed[region][2])))
    #             intersect = set_call_pos.intersection(set_region_pos)
    #             if len(intersect)>1:
    #                 results.append({"chr":call[0],"region":region,"start":int(call[1]),"end":int(call[2])})
    #     return results
