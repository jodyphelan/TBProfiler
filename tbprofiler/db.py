import csv
import json
import re
from collections import defaultdict
import sys
from datetime import datetime
from pathogenprofiler import run_cmd, cmd_out
from .utils import load_gff
import os
import shutil

chr_name = "Chromosome"

def fa2dict(filename):
    fa_dict = {}
    seq_name = ""
    for l in open(filename):
        line = l.rstrip()
        if line=="":continue
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


def write_bed(gene_dict,gene_info,outfile):
    lines = []
    for gene in gene_dict:
        if gene not in gene_info:
            sys.stderr.write("%s not found in the 'gene_info' dictionary... Exiting!" % gene)
            quit()
        lines.append([
            gene_info[gene].chrom,
            str(gene_info[gene].feature_start-200),
            str(gene_info[gene].feature_end+200),
            gene_info[gene].locus_tag,
            gene_info[gene].name,
            ",".join(gene_dict[gene])
        ])
    with open(outfile,"w") as O:
        for line in sorted(lines,key=lambda x: int(x[1])):
            O.write("%s\n" %"\t".join(line))

def load_gene_info(filename):
    gene_info = {}
    for l in open(filename):
        row = l.rstrip().split()
        strand = "-" if row[0][-1]=="c" else "+"
        gene_info[row[0]] = {"locus_tag":row[0],"gene":row[1],"start":int(row[2]),"end":int(row[3]),"gene_start":int(row[4]),"gene_end":int(row[5]),"strand":strand}
        gene_info[row[1]] = {"locus_tag":row[0],"gene":row[1],"start":int(row[2]),"end":int(row[3]),"gene_start":int(row[4]),"gene_end":int(row[5]),"strand":strand}
    return gene_info



def get_genome_position(gene_object,change):
    if change in ["frameshift","large_deletion","functional_gene"]:
        return None
    if "any_missense_codon" in change:
        codon = int(change.replace("any_missense_codon_",""))
        change = f"p.Xyz{codon}Xyz"


    g = gene_object
    r = re.search("p.[A-Za-z]+([0-9]+)",change)
    if r:
        codon = int(r.group(1))
        if g.strand=="+":
            p = g.start + codon*3-3
            return [p,p+1,p+2]
        else:
            p = g.start - codon*3 + 1
            return [p,p+1,p+2]
    r = re.search("c.(-[0-9]+)[ACGT]+>[ACGT]+",change)
    if r:
        pos = int(r.group(1))
        if g.strand=="+":
            p = g.start + pos
            return [p]
        else:
            p = g.start - pos
            return [p]
    r = re.search("n.([0-9]+)[ACGT]+>[ACGT]+",change)
    if r:
        pos = int(r.group(1))
        if g.strand=="+":
            p = g.start + pos -1
            return [p]
    r = re.search("c.([0-9]+)_([0-9]+)ins[A-Z]+",change)
    if r:
        pos = int(r.group(1))
        if g.strand=="+":
            p = g.start + pos -1
            return [p, p+1]
        else:
            p = g.start - pos 
            return [p, p+1]
    r = re.search("c.([\-0-9]+)_([\-0-9]+)del[A-Z]+",change)
    if r:
        pos1 = int(r.group(1))
        pos2 = int(r.group(2))
        if g.strand=="+":
            p1 = g.start + pos1 -1
            p2 = g.start + pos2 -1
            if pos1<0:
                p1+=1
                p2+=1
            return list(range(p1,p2+1))
        else:
            p1 = g.start - pos1 + 1
            p2 = g.start - pos2 + 1
            if pos1<0:
                p1+=1
                p2+=1
                quit(f"Don't know how to handle {change}")
            return list(range(p2,p1+1))
    r = re.search("c.([0-9]+)del[A-Z]+",change)
    if r:
        pos = int(r.group(1))
        if g.strand=="+":
            p = g.start + pos - 1
            return [p]
        else:
            p = g.start - pos + 1
            return [p]
    r = re.search("c.([0-9]+)dup[A-Z]+",change)
    if r:
        pos = int(r.group(1))
        if g.strand=="+":
            p = g.start + pos - 1
            return [p]
        else:
            p = g.start - pos + 1
            return [p]
    r = re.search("c.([0-9]+)_([0-9]+)dup[A-Z]+",change)
    if r:
        pos1 = int(r.group(1))
        pos2 = int(r.group(2))
        if g.strand=="+":
            p1 = g.start + pos1 - 1
            p2 = g.start + pos2 - 1
            return list(range(p1,p2+1))
        else:
            p = g.start - pos + 1
            quit(f"Don't know how to handle {change}")
            return [p]


    r = re.search("n.([0-9]+)([0-9]+)dup[A-Z]+",change)
    if r:
        pos1 = int(r.group(1))
        pos2 = int(r.group(2))
        if g.strand=="+":
            p1 = g.start + pos1 - 1
            p2 = g.start + pos2 - 1
            return list(range(p1,p2+1))
        else:
            p = g.start - pos + 1
            quit(f"Don't know how to handle {change}")
            return [p]
    quit(f"Don't know how to handle {change}")

def match_ref_chrom_names(source,target):
    source_fa = fa2dict(source)
    source_fa_size = {s:len(source_fa[s]) for s in source_fa}
    target_fa = fa2dict(target)
    target_fa_size = {s:len(target_fa[s]) for s in target_fa}
    conversion = {}
    for s in target_fa:
        tlen = target_fa_size[s]
        tmp = [x[0] for x in source_fa_size.items() if x[1]==tlen]
        if len(tmp)==1:
            conversion[s] = tmp[0]
    return conversion


def create_db(args):

    genome_file = "%s.fasta" % args.prefix
    gff_file = "%s.gff" % args.prefix
    if args.prefix:
        barcode_file = "%s.barcode.bed" % args.prefix
    bed_file = "%s.bed" % args.prefix
    json_file = "%s.dr.json" % args.prefix
    version_file = "%s.version.json" % args.prefix


    global chr_name

    if args.match_ref:
        chrom_conversion = match_ref_chrom_names(args.match_ref,"genome.fasta")
        shutil.copy(args.match_ref,genome_file)
    else:
        chrom_conversion = match_ref_chrom_names("genome.fasta","genome.fasta")
        shutil.copy("genome.fasta",genome_file)    
    
    with open(gff_file,"w") as O:
        for l in open("genome.gff"):
            if l[0]=="#":
                O.write(l)
            else:
                row = l.strip().split()
                if row[0] in chrom_conversion:
                    row[0] = chrom_conversion[row[0]]
                    O.write("\t".join(row)+"\n")        

                    
    genes = load_gff(gff_file)
    gene_name2gene_id = {g.name:g.locus_tag for g in genes.values()}
    gene_name2gene_id.update({g.locus_tag:g.locus_tag for g in genes.values()})
    db = {}
    locus_tag_to_drug_dict = defaultdict(set)
    for row in csv.DictReader(open(args.csv)):
        locus_tag = gene_name2gene_id[row["Gene"]]
        drug = row["Drug"].lower()
        mut = row["Mutation"]
        locus_tag_to_drug_dict[locus_tag].add(drug)
        if locus_tag not in db:
            db[locus_tag] = {}
        if mut not in db[locus_tag]:
            db[locus_tag][mut] = {"annotations":[]}

        tmp_annotation = {"type":"drug","drug":row["Drug"]}
        annotation_columns = set(row.keys()) - set(["Gene","Mutation","Drug"])
        for col in annotation_columns:
            if row[col]=="":continue
            tmp_annotation[col.lower()] = row[col]
        db[locus_tag][mut]["annotations"].append(tmp_annotation)
        db[locus_tag][mut]["genome_positions"] = get_genome_position(genes[locus_tag],mut)

    if args.other_annotations:
        for row in csv.DictReader(open(args.other_annotations)):
            locus_tag = gene_name2gene_id[row["Gene"]]
            mut = row["Mutation"]
            if locus_tag not in db:
                db[locus_tag] = {}
            if mut not in db[locus_tag]:
                db[locus_tag][mut] = {"annotations":[]}
            tmp_annotation = {"type":row["Type"]}
            annotation_columns = set(row.keys()) - set(["Gene","Mutation"])
            for col in annotation_columns:
                if row[col]=="":continue
                tmp_annotation[col.lower()] = row[col]
            db[locus_tag][mut]["annotations"].append(tmp_annotation)
            db[locus_tag][mut]["genome_positions"] = get_genome_position(genes[locus_tag],mut)

    if args.watchlist:
        for row in csv.DictReader(open(args.watchlist)):
            locus_tag = gene_name2gene_id[row["Gene"]]
            drug = row["Drug"].lower()
            locus_tag_to_drug_dict[locus_tag].add(drug)

    

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

    variables = json.load(open("variables.json"))    
    variables["chromosome_conversion"] = {"source":list(chrom_conversion.keys()),"target":list(chrom_conversion.values())}
    json.dump(variables,open(args.prefix+".variables.json","w"))
    if os.path.isfile("barcode.bed"):
        with open(barcode_file,"w") as O:
            for l in open("barcode.bed"):
                row = l.strip().split()
                row[0] = chrom_conversion[row[0]]
                O.write("\t".join(row)+"\n")
    write_bed(locus_tag_to_drug_dict,genes,bed_file)
    json.dump(db,open(json_file,"w"))
