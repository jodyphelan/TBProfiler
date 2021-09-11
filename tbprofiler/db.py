import csv
import json
import re
from collections import defaultdict
import sys
from datetime import datetime
from pathogenprofiler import run_cmd, cmd_out

chr_name = "Chromosome"

def fa2dict(filename):
    fa_dict = {}
    seq_name = ""
    for l in open(filename):
        line = l.rstrip()
        if line[0] == ">":
            seq_name = line[1:].split()[0]
            fa_dict[seq_name] = []
        else:
            fa_dict[seq_name].append(line)
    result = {}
    for seq in fa_dict:
        result[seq] = "".join(fa_dict[seq])
    return result

def revcom(s):
    """Return reverse complement of a sequence"""
    def complement(s):
                    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
                    letters = list(s)
                    letters = [basecomplement[base] for base in letters]
                    return ''.join(letters)
    return complement(s[::-1])

def write_gene_pos(infile,genes,outfile):
    with open(outfile, "w") as OUT:
        for l in open(infile):
            row = l.strip().split()
            rv,_,chr_start,chr_end,gene_start,gene_end = [row[0],row[1]]+[int(row[i]) for i in range(2,6)]
            if rv in genes:
                y = 0
                for i, chr_pos in enumerate(range(chr_start, chr_end+1)):
                    x = 1 if gene_start< gene_end else -1
                    if gene_start+(x*i) == 0:
                        y = 1 if gene_start< gene_end else -1
                    OUT.write("%s\t%s\t%s\t%s\n" % (chr_name,chr_pos,rv,gene_start+(x*i)+y))


def write_bed(gene_dict,gene_info,outfile,chr_name):
    O = open(outfile,"w")
    lines = []
    for gene in gene_dict:
        if gene not in gene_info:
            sys.stderr.write("%s not found in the 'gene_info' dictionary... Exiting!" % gene)
            quit()
        lines.append([chr_name,int(gene_info[gene]["start"]),int(gene_info[gene]["end"]),gene_info[gene]["locus_tag"],gene_info[gene]["gene"],",".join(gene_dict[gene])])
    for line in sorted(lines,key=lambda x: x[1]):
        line[1] = str(line[1])
        line[2] = str(line[2])
        O.write("%s\n" %"\t".join(line))
    O.close()

def load_gene_info(filename):
    gene_info = {}
    for l in open(filename):
        row = l.rstrip().split()
        strand = "-" if row[0][-1]=="c" else "+"
        gene_info[row[0]] = {"locus_tag":row[0],"gene":row[1],"start":int(row[2]),"end":int(row[3]),"gene_start":int(row[4]),"gene_end":int(row[5]),"strand":strand}
        gene_info[row[1]] = {"locus_tag":row[0],"gene":row[1],"start":int(row[2]),"end":int(row[3]),"gene_start":int(row[4]),"gene_end":int(row[5]),"strand":strand}
    return gene_info

def create_db(args):
    global chr_name
    chr_name = args.seqname
    fasta_dict = fa2dict("genome.fasta")
    gene_info = load_gene_info("genes.txt")
    db = {}
    locus_tag_to_drug_dict = defaultdict(set)
    confidence = {}
    for row in csv.DictReader(open(args.confidence)):
        confidence[(row["gene"],row["mutation"],row["drug"])] = row["confidence"]
    for row in csv.DictReader(open(args.csv)):
        locus_tag = gene_info[row["Gene"]]["locus_tag"]
        drug = row["Drug"].lower()
        mut = row["Mutation"]
        locus_tag_to_drug_dict[locus_tag].add(drug)
        if locus_tag not in db:
            db[locus_tag] = {}
        if mut not in db[locus_tag]:
            db[locus_tag][mut] = {"annotations":[]}

        tmp_annotation = {"type":"drug","drug":row["Drug"]}
        tmp_annotation["confidence"] = confidence.get((locus_tag,row["Mutation"],drug),"indeterminate")
        annotation_columns = set(row.keys()) - set(["Gene","Mutation","Drug"])
        for col in annotation_columns:
            if row[col]=="":continue
            tmp_annotation[col.lower()] = row[col]
        db[locus_tag][mut]["annotations"].append(tmp_annotation)


    for row in csv.DictReader(open(args.watchlist)):
        locus_tag = gene_info[row["Gene"]]["locus_tag"]
        drug = row["Drug"].lower()
        locus_tag_to_drug_dict[locus_tag].add(drug)

    genome_file = "%s.fasta" % args.prefix
    gff_file = "%s.gff" % args.prefix
    ann_file = "%s.ann.txt" % args.prefix
    barcode_file = "%s.barcode.bed" % args.prefix
    bed_file = "%s.bed" % args.prefix
    json_file = "%s.dr.json" % args.prefix
    version_file = "%s.version.json" % args.prefix

    version = {"name":args.prefix}
    if not args.custom:
        for l in cmd_out("git log | head -4"):
            row = l.strip().split()
            if row == []: continue
            version[row[0].replace(":","")] = " ".join(row[1:])
        version["commit"] = version["commit"][:7]
    else:
        version["Date"] = str(datetime.now()) if not args.db_date else args.db_date
        version["name"] = args.db_name if args.db_name else "NA"
        version["commit"] = args.db_commit if args.db_name else "NA"
        version["Author"] = args.db_author if args.db_author else "NA"

    json.dump(version,open(version_file,"w"))
    open(genome_file,"w").write(">%s\n%s\n" % (chr_name,fasta_dict["Chromosome"]))
    run_cmd("sed 's/Chromosome/%s/g' genome.gff > %s" % (chr_name,gff_file))
    run_cmd("sed 's/Chromosome/%s/g' barcode.bed > %s" % (chr_name,barcode_file))
    write_gene_pos("genes.txt",list(locus_tag_to_drug_dict.keys()),ann_file)
    write_bed(locus_tag_to_drug_dict,gene_info,bed_file,chr_name)
    json.dump(db,open(json_file,"w"))
