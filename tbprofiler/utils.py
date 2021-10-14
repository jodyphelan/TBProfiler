import json
import re
from collections import defaultdict
import sys
import pathogenprofiler as pp
import os

def get_conf_dict_with_path(library_path):
    files = {"gff":".gff","ref":".fasta","barcode":".barcode.bed","bed":".bed","json_db":".dr.json","version":".version.json"}
    conf = {}
    for key in files:
        sys.stderr.write("Using %s file: %s\n" % (key,library_path+files[key]))
        if os.path.isfile(library_path+files[key]):
            conf[key] = pp.filecheck(library_path+files[key])
#     test = json.load(open(conf["json_db"]))["Rv1908c"]["p.Ser315Thr"]
#     if "annotation" not in test and "drugs" in test:
#         quit("""\n
# ################################# ERROR #######################################

# The database has different format than expected. Since tb-profiler v2.4 the
# database is parsed using tb-profiler code. Please run the following code to get
# the latest version of the database or load your own:

# tb-profiler update_tbdb

# or

# tb-profiler load_library /path/to/custom_library

# ###############################################################################
# """)


    return conf

def get_conf_dict(library_prefix):
    library_prefix = "%s/share/tbprofiler/%s" % (sys.base_prefix,library_prefix)
    return get_conf_dict_with_path(library_prefix)

def get_lt2drugs(bed_file):
    lt2drugs = {}
    for l in open(bed_file):
        row = l.strip().split()
        lt2drugs[row[3]] = row[5].split(",")
    return lt2drugs

def get_gene2drugs(bed_file):
    lt2drugs = {}
    for l in open(bed_file):
        row = l.strip().split()
        lt2drugs[row[4]] = row[5].split(",")
    return lt2drugs

def get_drugs2lt(bed_file):
    tmp = get_lt2drugs(bed_file)
    results = defaultdict(list)
    for gene in tmp:
        for drug in tmp[gene]:
            results[drug].append(gene)
    return dict(results)

def get_drugs2gene(bed_file):
    tmp = get_gene2drugs(bed_file)
    results = defaultdict(list)
    for gene in tmp:
        for drug in tmp[gene]:
            results[drug].append(gene)
    return dict(results)

class gene_class:
    def __init__(self,name,locus_tag,strand,chrom,start,end,length):
        self.name = name
        self.locus_tag = locus_tag
        self.strand = strand
        self.chrom = chrom
        self.feature_start = start
        self.feature_end = end
        self.start = self.feature_start if strand=="+" else self.feature_end
        self.end = self.feature_end if strand=="+" else self.feature_start
        self.length = length

def load_gff(gff):
    genes = {}
    for l in open(gff):
        if l[0]=="#": continue
        fields = l.rstrip().split()
        if fields[2] not in ["gene","rRNA_gene","ncRNA_gene"]: continue
        strand = fields[6]
        chrom = fields[0]
        p1 = int(fields[3])
        p2 = int(fields[4])
        gene_length = p2-p1+1
        re_obj = re.search("Name=([a-zA-Z0-9\.\-\_\(\)]+)",l)
        gene_name = re_obj.group(1) if re_obj else "NA"
        re_obj = re.search("gene_id=([a-zA-Z0-9\.\-\_]+)",l)
        locus_tag = re_obj.group(1) if re_obj else "NA"
        start = p1
        end =  p2
        tmp = gene_class(gene_name,locus_tag,strand,chrom,start,end,gene_length)
        genes[locus_tag] = tmp
    return genes

def get_genome_positions_from_json_db(json_file):
    genome_positions = defaultdict(set)
    db = json.load(open(json_file))
    for gene in db:
        for var in db[gene]:
            drugs = tuple([x["drug"] for x in db[gene][var]["annotations"]])
            if db[gene][var]["genome_positions"]:
                for pos in db[gene][var]["genome_positions"]:
                    genome_positions[pos].add((gene,var,drugs))

    return genome_positions


def rv2genes(bed_file):
    #Chromosome      759310  763325  Rv0667  rpoB    rifampicin
    rv2gene = {}
    for l in open(bed_file):
        row = l.strip().split()
        rv2gene[row[3]] = row[4]
    return rv2gene
