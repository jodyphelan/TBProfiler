import csv
import json
import re
from collections import defaultdict
import sys
from datetime import datetime
from pathogenprofiler import run_cmd, cmd_out, errlog, unlist, debug
from .utils import load_gff
import os
import shutil
from uuid import uuid4
import pathogenprofiler as pp

chr_name = "Chromosome"

def generate_kmer_database(kmer_file,outfile):
    from itertools import combinations, product

    def generate(s, d=1):
        N = len(s)
        letters = 'ACGT'
        pool = list(s)

        for indices in combinations(range(N), d):
            for replacements in product(letters, repeat=d):
                skip = False
                for i, a in zip(indices, replacements):
                    if pool[i] == a: skip = True
                if skip: continue

                keys = dict(zip(indices, replacements))
                yield ''.join([pool[i] if i not in indices else keys[i] 
                            for i in range(N)])

    with open(outfile,"w") as O:
        for l in open(kmer_file):
            row = l.strip().split()
            kmers = [row[0]] + list(generate(row[0]))
            O.write("%s\t%s\t%s\n" % (row[1],len(kmers),"\t".join(kmers)))



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

def extract_genome_positions(db,gene):
    pos = []
    for mut in db[gene]:
        if any([a["type"]=="drug" for a in db[gene][mut]["annotations"]]):
            if mut in ["functional_gene","frameshift","large_deletion","transcript_ablation"]: continue
            pos.extend(db[gene][mut]["genome_positions"])
    return list(set(pos))

def write_bed(db,gene_dict,gene_info,outfile,padding=200):
    lines = []
    for gene in gene_dict:
        if gene not in gene_info:
            errlog("%s not found in the 'gene_info' dictionary... Exiting!" % gene)
            quit()
        if gene_info[gene].locus_tag in db:
            genome_positions = extract_genome_positions(db,gene_info[gene].locus_tag)
            if len(genome_positions)>0 and (gene_info[gene].feature_start > min(genome_positions)):
                genome_start = min(genome_positions) - padding
            else:
                genome_start = gene_info[gene].feature_start - padding
            
            if len(genome_positions)>0 and (gene_info[gene].feature_end < max(genome_positions)):
                genome_end = max(genome_positions) + padding
            else:
                genome_end = gene_info[gene].feature_end + padding
        else:
            genome_start = gene_info[gene].feature_start - padding
            genome_end = gene_info[gene].feature_end + padding

        lines.append([
            gene_info[gene].chrom,
            str(genome_start),
            str(genome_end),
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

def get_ann(variants,snpEffDB):
    uuid = str(uuid4()) #"463545ef-71fc-449b-8f4e-9c907ee6fbf5"
    with open(uuid,"w") as O:
        O.write('##fileformat=VCFv4.2\n')
        O.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        O.write('##contig=<ID=Chromosome,length=4411532>\n')
        O.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttest\n')
        for var in variants.values():
            O.write("Chromosome\t%(pos)s\t.\t%(ref)s\t%(alt)s\t255\t.\t.\tGT\t1\n" % var) 
    results = {}
    keys = list(variants.keys())
    vals = list(variants.values())
    i = 0
    for l in cmd_out(f"snpEff ann {snpEffDB} {uuid}"):
        if l[0]=="#": continue
        row = l.strip().split()
        for ann in row[7].split(","):
            a = ann.split("|")
            if len(a)!=16:continue

            if vals[i]["gene"] in [a[3],a[4]]:
                results[keys[i]] = a[9] if vals[i]["type"]=="nucleotide" else a[10]
        i+=1
    os.remove(uuid)
    return results


def get_snpeff_formated_mutation_list(csv_file,ref,gff,snpEffDB):
    genes = load_gff(gff,aslist=True)
    refseq = fa2dict(ref)

    mutations  =  {}
    converted_mutations = {}
    for row in csv.DictReader(open(csv_file)):
        gene = [g for g in genes if g.name==row["Gene"] or g.locus_tag==row["Gene"]][0]
        r = re.search("n.([0-9]+)([ACGT]+)>([ACGT]+)",row["Mutation"])
        if r:
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
            
        r = re.search("n.([0-9]+)([ACGT]+)>([ACGT]+)",row["Mutation"])
        if r:
            converted_mutations[(row["Gene"],row["Mutation"])] = f"n.{r.group(1)}{r.group(2).upper()}>{r.group(3).upper()}"
        r = re.search("p\..+",row["Mutation"])
        if r:
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
        r = re.search("c.-[0-9]+[ACGT]>[ACGT]",row["Mutation"])
        if r:
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]

        r = re.search("c.[0-9]+dup[ACGT]+",row["Mutation"])
        if r:
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
        r = re.search("c.[0-9]+_[0-9]+dup[ACGT]+",row["Mutation"])
        if r:
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
        



        r = re.search("c.([0-9]+)del",row["Mutation"])
        if r:
            # "ethA" "c.341del"
            del_start = int(r.group(1))
            del_end = int(r.group(1))
            if gene.strand == "+":
                # rpoB "c.1282_1290del"
                genome_start = gene.start + del_start - 2
                genome_end = gene.start + del_end 
            else:
                # "ethA" "c.1057_1059del"
                genome_start = gene.start - del_end
                genome_end = gene.start - del_start + 2
            ref = refseq["Chromosome"][genome_start-1:genome_end-1]
            alt = ref[0]
            mutations[(row["Gene"],row["Mutation"])] = {"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}

        r = re.search("c.([0-9]+)_([0-9]+)del",row["Mutation"])
        if r:
            del_start = int(r.group(1))
            del_end = int(r.group(2))
            if gene.strand == "+":
                # rpoB "c.1282_1290del"
                genome_start = gene.start + del_start - 2
                genome_end = gene.start + del_end 
            else:
                # "ethA" "c.1057_1059del"
                genome_start = gene.start - del_end
                genome_end = gene.start - del_start + 2
            ref = refseq["Chromosome"][genome_start-1:genome_end-1]
            alt = ref[0]
            mutations[(row["Gene"],row["Mutation"])] = {"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}

        r = re.search("c.-([0-9]+)del",row["Mutation"])
        if r:
            del_start = int(r.group(1))
            del_end = int(r.group(1))
            if gene.strand == "+":
               # "embA" "c.-29_-28del"
                genome_start = gene.start - del_start - 1
                genome_end = gene.start - del_end + 1
            else:
                # "alr" "c.-283_-280delCAAT"
                genome_start = gene.start + del_end - 1
                genome_end = gene.start + del_start + 1

            ref = refseq["Chromosome"][genome_start-1:genome_end-1]
            alt = ref[0]
            mutations[(row["Gene"],row["Mutation"])] = {"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}

        
        r = re.search("c.(-[0-9]+)_(-[0-9]+)del",row["Mutation"])
        if r:
            del_start = int(r.group(1))
            del_end = int(r.group(2))
            if gene.strand == "+":
               # "embA" "c.-29_-28del"
                genome_start = gene.start + del_start - 1
                genome_end = gene.start + del_end + 1
            else:
                # "alr" "c.-283_-280delCAAT"
                genome_start = gene.start - del_end - 1
                genome_end = gene.start - del_start + 1
            ref = refseq["Chromosome"][genome_start-1:genome_end-1]
            alt = ref[0]
            mutations[(row["Gene"],row["Mutation"])] = {"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}

        
        r = re.search("c.(-[0-9]+)_([0-9]+)del",row["Mutation"])
        if r:
            del_start = int(r.group(1))
            del_end = int(r.group(2))
            if gene.strand == "+":
               # "ethA" "c.-1058_968del"
                genome_start = gene.start + del_start -1
                genome_end = gene.start + del_end 
                # quit("Need to define!")

            else:
               # "ethA" "c.-1058_968del"
                genome_start = gene.start - del_end 
                genome_end = gene.start - del_start + 1

            ref = refseq["Chromosome"][genome_start-1:genome_end-1]
            alt = ref[0]
            mutations[(row["Gene"],row["Mutation"])] = {"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}


        r = re.search("c.([0-9]+)_([0-9]+)ins([ACGT]+)", row["Mutation"])
        if r:
            ins_start = int(r.group(1))
            ins_end = int(r.group(2))
            ins_seq = r.group(3)
            if gene.strand == "+":
                # "rpoB" "c.1296_1297insTTC"
                genome_start = gene.start + ins_start - 1 
                genome_end = gene.start + ins_end - 1
            else:
                # "pncA" "c.521_522insT"
                ins_seq = pp.revcom(ins_seq)
                genome_start = gene.start - ins_start 
                genome_end = gene.start - ins_end + 2

            ref = refseq["Chromosome"][genome_start-1:genome_end-1]
            alt = ref + ins_seq
            mutations[(row["Gene"],row["Mutation"])] = {"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}
        

        r = re.search("c.(-[0-9]+)_(-[0-9]+)ins([ACGT]+)",row["Mutation"])
        if r:
            del_start = int(r.group(1))
            del_end = int(r.group(2))
            ins_seq = r.group(3)
            if gene.strand == "+":
               # "rrs" "c.-29_-28insATAC"
                genome_start = gene.start + del_start 
                genome_end = gene.start + del_end 
            else:
                # "alr" "c.-283_-280delCAAT"
                ins_seq = pp.revcom(ins_seq)
                genome_start = gene.start - del_end 
                genome_end = gene.start - del_start 
            ref = refseq["Chromosome"][genome_start-1:genome_end-1]
            alt = ref + ins_seq
            mutations[(row["Gene"],row["Mutation"])] = {"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}


        if row["Mutation"] == "frameshift":
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
        if row["Mutation"] == "large_deletion":
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
        if row["Mutation"] == "transcript_ablation":
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
        if row["Mutation"] == "functional_gene":
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
        if row["Mutation"][:19] == "any_missense_codon_":
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
        
        if (row["Gene"],row["Mutation"]) not in converted_mutations and (row["Gene"],row["Mutation"]) not in mutations:
            quit(f"Don't know how to handle this mutation: {row['Gene']} {row['Mutation']}\n")
            
    print("Converting %s mutations" % len(mutations))
    if len(mutations)>0:
        mutation_conversion = get_ann(mutations,snpEffDB)
        for key in mutation_conversion:
            converted_mutations[key] = mutation_conversion[key]
    return converted_mutations
    



def get_genome_position(gene_object,change):
    if change in ["frameshift","large_deletion","functional_gene","transcript_ablation"]:
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
    r = re.search("[nc].([\-\*0-9]+)_([\-\*0-9]+)ins[A-Z]+",change)
    if r:
        if "*" in r.group(1):
            if g.strand=="+":
                pos = g.feature_end - g.feature_start + int(r.group(1).replace("*","")) + 1
            else:
                pos = g.feature_end - g.feature_start - int(r.group(1).replace("*","")) - 1
        else:
            pos = int(r.group(1))
        if g.strand=="+":
            p = g.start + pos -1
            return [p, p+1]
        else:
            p = g.start - pos 
            return [p, p+1]

    r = re.search("[nc].([\-\*0-9]+)_([\-\*0-9]+)del[A-Z]*",change)
    if r:
        if "*" in r.group(1):
            if g.strand=="+":
                pos1 = g.feature_end - g.feature_start + int(r.group(1).replace("*",""))
            else:
                pos1 = g.feature_end - g.feature_start - int(r.group(1).replace("*",""))
        else:
            pos1 = int(r.group(1))
        if "*" in r.group(2):
            if g.strand=="+":
                pos2 = g.feature_end - g.feature_start + int(r.group(2).replace("*",""))
            else:
                pos2 = g.feature_end - g.feature_start - int(r.group(2).replace("*",""))
        else:
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
                p1-=1
            if pos2<0:
                p2-=1
            return list(range(p2,p1+1))

    r = re.search("[nc].([\-\*0-9]+)del[A-Z]+",change)
    if r:
        if "*" in r.group(1):
            if g.strand=="+":
                pos = g.feature_end - g.feature_start + int(r.group(1).replace("*","")) + 1
            else:
                pos = g.feature_end - g.feature_start - int(r.group(1).replace("*","")) - 1
        else:
            pos = int(r.group(1))
        if g.strand=="+":
            p = g.start + pos - 1
            return [p]
        else:
            p = g.start - pos + 1
            return [p]

    r = re.search("[nc].([\-0-9]+)dup[A-Z]+",change)
    if r:
        pos = int(r.group(1))
        if g.strand=="+":
            p = g.start + pos - 1
            return [p]
        else:
            p = g.start - pos + 1
            return [p]
    
    r = re.search("[nc].([\-0-9]+)_([\-0-9]+)dup[A-Z]+",change)
    if r:
        pos1 = int(r.group(1))
        pos2 = int(r.group(2))
        if g.strand=="+":
            p1 = g.start + pos1 - 1
            p2 = g.start + pos2 - 1
            return list(range(p1,p2+1))
        else:
            p1 = g.start - pos1 + 1
            p2 = g.start - pos2 + 1
            return list(range(p2,p1+1))

    

    # r = re.search("n.([0-9]+)([0-9]+)dup[A-Z]+",change)
    # if r:
    #     pos1 = int(r.group(1))
    #     pos2 = int(r.group(2))
    #     if g.strand=="+":
    #         p1 = g.start + pos1 - 1
    #         p2 = g.start + pos2 - 1
    #         return list(range(p1,p2+1))
    #     else:
    #         p = g.start - pos + 1
    #         quit(f"Don't know how to handle {change}")
    #         return [p]
    quit(f"Don't know how to handle {str(vars(g))} {change}")

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
    mutation_lookup = get_snpeff_formated_mutation_list(args.csv,"genome.fasta","genome.gff",json.load(open("variables.json"))["snpEff_db"])
    with open(args.prefix+".conversion.log","w") as L:
        for row in csv.DictReader(open(args.csv)):
            locus_tag = gene_name2gene_id[row["Gene"]]
            drug = row["Drug"].lower()
            mut = mutation_lookup[(row["Gene"],row["Mutation"])]
            if args.include_original_mutation:
                row["original_mutation"] = row["Mutation"]
            if mut!=row["Mutation"]:
                L.write(f"Converted {row['Gene']} {row['Mutation']} to {mut}\n")
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
            mutation_lookup = get_snpeff_formated_mutation_list(args.other_annotations,"genome.fasta","genome.gff",json.load(open("variables.json"))["snpEff_db"])
            for row in csv.DictReader(open(args.other_annotations)):
                locus_tag = gene_name2gene_id[row["Gene"]]
                mut = mutation_lookup[(row["Gene"],row["Mutation"])]
                if mut!=row["Mutation"]:
                    L.write(f"Converted {row['Gene']} {row['Mutation']} to {mut}\n")
                if locus_tag not in db:
                    db[locus_tag] = {}
                if mut not in db[locus_tag]:
                    db[locus_tag][mut] = {"annotations":[]}
                tmp_annotation = {"type":row["Type"]}
                if args.include_original_mutation:
                    tmp_annotation["original_mutation"] = row["Mutation"]


                for x in row["Info"].split(";"):
                    key,val = x.split("=")
                    tmp_annotation[key.lower()] = val
                    if key=="drug":
                        locus_tag_to_drug_dict[locus_tag].add(val)
                db[locus_tag][mut]["annotations"].append(tmp_annotation)
                db[locus_tag][mut]["genome_positions"] = get_genome_position(genes[locus_tag],mut)

        if args.watchlist:
            for row in csv.DictReader(open(args.watchlist)):
                locus_tag = gene_name2gene_id[row["Gene"]]
                for d in row["Drug"].split(","):
                    drug = d.lower()
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
        variables["chromosome_conversion"] = {"target":list(chrom_conversion.keys()),"source":list(chrom_conversion.values())}
        json.dump(variables,open(args.prefix+".variables.json","w"))
        if os.path.isfile("barcode.bed"):
            with open(barcode_file,"w") as O:
                for l in open("barcode.bed"):
                    row = l.strip().split("\t")
                    row[0] = chrom_conversion[row[0]]
                    O.write("\t".join(row)+"\n")
        write_bed(db,locus_tag_to_drug_dict,genes,bed_file)
        json.dump(db,open(json_file,"w"))
        if args.spoligotypes:
            spoligotype_file = args.prefix+".spoligotype_spacers.txt"
            generate_kmer_database("spoligotype_spacers.txt", spoligotype_file)