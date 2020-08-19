from .utils import *
from .fasta import fasta
from collections import defaultdict
import re
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
    def csq(self,ref_file,gff_file):
        add_arguments_to_self(self,locals())
        self.vcf_csq_file = self.prefix+".csq.vcf.gz"
        run_cmd("bcftools csq -p m -f %(ref_file)s -g %(gff_file)s %(filename)s -Oz -o %(vcf_csq_file)s" % vars(self))
        return vcf(self.vcf_csq_file,self.prefix)

    def add_annotations(self,ref_file,bam_file):
        add_arguments_to_self(self,locals())
        self.new_file = self.prefix + ".ann.vcf.gz"

        run_cmd("gatk VariantAnnotator -R %(ref_file)s -I %(bam_file)s -V %(filename)s -O %(new_file)s -A MappingQualityRankSumTest -A ReadPosRankSumTest -A QualByDepth -A BaseQualityRankSumTest -A TandemRepeat -A StrandOddsRatio -OVI false" % vars(self))
        return vcf(self.new_file,self.prefix)

    def load_csq(self, ann_file = None, allow_intergenic=False,debug = False, extract_ann=False):
        ann = defaultdict(dict)
        gene_set = set()
        if ann_file:
            for l in open(ann_file):
                #chrom pos gene gene/codon_pos
                row = l.rstrip().split()
                ann[row[0]][int(row[1])] = (row[2],row[3])
                gene_set.add(row[2])

        nuc_variants = self.load_variants()
        variants = {s:[] for s in self.samples}
        annotations_types = self.get_gatk_annotations() if extract_ann else []
        annotation_extract_string = "\t%s" % "\\t".join(["%%%s" % x for x in annotations_types]) if extract_ann else ""
        csq_start_column = 4 + len(annotations_types)
        for line in cmd_out("bcftools query -u -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT%s[\\t%%SAMPLE\\t%%TBCSQ\\t%%TGT\\t%%AD]\\n' %s" % (annotation_extract_string,self.filename)):

            if debug:
                sys.stderr.write(line+"\n")
            row = line.split()
            chrom = row[0]
            pos = int(row[1])
            ref = row[2]
            alts = row[3].split(",")
            alleles = [ref]+alts
            annotations = {}
            for i,j in enumerate(range(4,csq_start_column)):
                annotations[annotations_types[i]] = row[j]
            if chrom in ann and pos in ann[chrom]:
                ann_pos = int(ann[chrom][pos][1])
                ann_gene = ann[chrom][pos][0]
            else:
                ann_pos = None
                ann_gene = None
            if len(row)==csq_start_column or "coding_sequence" in line:
                for alt in alts:
                    if chrom in ann and pos in ann[chrom]:
                        cng = "%s%s>%s" % (ann_pos,ref,alt)
                        nucleotide_change = "%s%s>%s" % (pos,ref,alt)
                        for sample in self.samples:
                            if sample in nuc_variants[chrom][pos] and alt in nuc_variants[chrom][pos][sample]:
                                variants[sample].append({"sample":sample,"gene_id":ann_gene,"chr":chrom,"genome_pos":pos,"type":"non_coding","change":cng,"nucleotide_change":nucleotide_change,"freq":nuc_variants[chrom][pos][sample][alt], "variant_annotations":annotations})
                    elif allow_intergenic:
                        cng = "%s%s>%s" % (pos,ref,alt)
                        nucleotide_change = "%s%s>%s" % (pos,ref,alt)
                        for sample in self.samples:
                            variants[sample].append({"sample":sample,"gene_id":"NA","chr":chrom,"genome_pos":pos,"type":"intergenic","change":cng,"nucleotide_change":nucleotide_change,"freq":nuc_variants[chrom][pos][sample][alt], "variant_annotations":annotations})
                    else:
                        log(line)
                        log("ERROR in loading alts",True)
                continue

            for i in range(csq_start_column,len(row)-4,5):
                # if pos==1472723:
                    # import pdb; pdb.set_trace()
                # if row[i+1] == ".": continue
                if row[i+1][0] == "@": continue

                sample = row[i]
                infos = [x.split("|") for x in row[i+1].split(",") if x!="."] + [x.split("|") for x in row[i+2].split(",") if x!="."]
                if debug:
                    sys.stderr.write(str(infos))
                    sys.stderr.write("\n")
                info = None
                for x in infos:
                    if x[1] in gene_set or x[2] in gene_set:
                        info = x
                        break
                if not ann_file:
                    info = infos[0]
                call1,call2 = row[i+3].split("/") if "/" in row[i+3] else row[i+3].split("|")
                ad = [int(x) if x!="." else 0 for x in row[i+4].split(",")]  if row[i+4]!="." else [0,100]

                adr = {alleles[i]:d/sum(ad) for i,d in enumerate(ad)}
                if info[-1] == "pseudogene": continue
                gene_name = info[1]
                gene_id = info[2] if info[2]!="" else gene_name
                if info[0] == "intron":continue
                if info[0] == "coding_sequence":
                    cng = "%s%s>%s" % (ann_pos,call1,call2)
                    variants[sample].append({"sample":sample,"gene_id":ann_gene,"chr":chrom,"genome_pos":pos,"type":"non_coding","change":cng,"freq":adr[call2], "nucleotide_change":cng, "variant_annotations":annotations})
                elif  "missense" in info[0] or "start_lost" in info[0] or "stop_gained" in info[0]:
                    variants[sample].append({"sample":sample,"gene_id":gene_id,"gene_name":gene_name,"chr":chrom,"genome_pos":pos,"type":info[0],"change":info[5],"freq":adr[call2], "nucleotide_change":info[6], "variant_annotations":annotations})
                elif "frame" in info[0] or "stop_lost" in info[0]:
                    if len(info)<6:
                        if chrom in ann and pos in ann[chrom]:
                            change = "%s%s>%s" % (pos,ref,call2)
                            variants[sample].append({"sample":sample,"gene_id":gene_id,"gene_name":gene_name,"chr":chrom,"genome_pos":pos,"type":info[0],"change":change,"freq":adr[call2],"nucleotide_change":change, "variant_annotations":annotations})
                    else:
                        variants[sample].append({"sample":sample,"gene_id":gene_id,"gene_name":gene_name,"chr":chrom,"genome_pos":pos,"type":info[0],"change":info[6],"freq":adr[call2],"nucleotide_change":info[6], "variant_annotations":annotations})
                elif "synonymous" in info[0] or info[0] == "stop_retained":
                    change_num,ref_nuc,alt_nuc =  parse_mutation(info[6])
                    change = "%s%s>%s" % (ann_pos,ref_nuc,alt_nuc) if ann_pos else "%s%s>%s" % (pos,ref_nuc,alt_nuc)
                    variants[sample].append({"sample":sample,"gene_id":gene_id,"gene_name":gene_name,"chr":chrom,"genome_pos":pos,"type":info[0],"change":change,"freq":adr[call2],"nucleotide_change":info[6], "variant_annotations":annotations})
                elif info[0] == "non_coding" or info[0] == "splice_region" or info[0] == "3_prime_utr":
                    if chrom in ann and pos in ann[chrom]:
                        gene = ann[chrom][pos][0]
                        gene_pos = ann[chrom][pos][1]
                        change = "%s%s>%s" % (gene_pos,ref,call2)
                        variants[sample].append({"sample":sample,"gene_id":gene,"gene_name":gene_name,"chr":chrom,"genome_pos":pos,"type":info[0],"change":change,"freq":adr[call2],"nucleotide_change":change, "variant_annotations":annotations})
                else:
                    log(line)
                    log(info[0]+"\n")
                    log("Unknown variant type...Exiting!\n",True)
        return variants

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
                call1,call2 = calls.split("/") if "/" in calls else calls.split("|")
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
    def get_robust_calls(self):
        results = []
        for l in cmd_out("bcftools query -f '%%CHROM\\t%%POS\\t[%%END\\t%%GT\\t%%DR\\t%%DV\\t%%RR\\t%%RV]\\n' %(filename)s" % vars(self)):
            row = l.split()
            if row[3]!="1/1":continue
            if int(row[2])-int(row[1])>100000: continue
            results.append(row)
        return results
    def overlap_bed(self,bed_file):
        results = []
        bed = load_bed(bed_file,[1,2,3,4,5],4)
        calls = self.get_robust_calls()
        for call in calls:
            set_call_pos = set(range(int(call[1]),int(call[2])))
            for region in bed:
                if bed[region][0]!=call[0]: continue
                set_region_pos = set(range(int(bed[region][1]),int(bed[region][2])))
                intersect = set_call_pos.intersection(set_region_pos)
                if len(intersect)>1:
                    results.append({"chr":call[0],"region":region,"start":min(intersect),"end":max(intersect)})
        return results
