import json
import re
from collections import defaultdict


def get_lt2drugs(bed_file):
    lt2drugs = {}
    for l in open(bed_file):
        row = l.strip().split()
        lt2drugs[row[3]] = row[5].split(",")
    return lt2drugs


def get_genome_positions_from_json_db(json_file,ann_file):
    codon_ann = defaultdict(set)
    gene_ann = {}
    for l in open(ann_file):
        #Chromosome      5002    Rv0005  -238
        chrom,genome_pos,gene,gene_pos = l.strip().split()
        genome_pos = int(genome_pos)
        gene_pos = int(gene_pos)
        # This is a bit dangerous as it assumes all coding genes are labled with "Rv"
        coding = True if gene[:2]=="Rv" else False
        if coding:
            if gene_pos>=0:
                codon = ((gene_pos-1)//3) + 1
                codon_ann[(gene,codon)].add(genome_pos)


        gene_ann[(gene,gene_pos)] = genome_pos


    genome_positions = defaultdict(set)
    db = json.load(open(json_file))
    for gene in db:
        for _var in db[gene]:
            var = db[gene][_var]["hgvs_mutation"]
            if var[0]=="p":
                #p.Thr40Ile
                re_obj = re.match("p.[A-Za-z]+([0-9]+)[A-Za-z\*]",var)
                codon_pos = int(re_obj.group(1))
                for pos in codon_ann[(gene,codon_pos)]:
                    genome_positions[pos].add((gene,var))
            elif var[0]=="c":
                if "ins" in var:
                    #c.192_193insG
                    re_obj = re.match("c.([0-9]+)_([0-9]+)ins[A-Za-z]+",var)
                    gene_pos = int(re_obj.group(1))
                    genome_positions[gene_ann[(gene,gene_pos)]].add((gene,var))

                elif "del" in var:
                    if "_" in var:
                        #c.785_787del
                        re_obj = re.match("c.(\-*[0-9]+)_(\-*[0-9]+)del",var)
                        gene_pos_from = int(re_obj.group(1))
                        gene_pos_to = int(re_obj.group(2))
                    else:
                        #c.884del
                        re_obj = re.match("c.(\-*[0-9]+)del",var)
                        gene_pos_from = int(re_obj.group(1))
                        gene_pos_to = gene_pos_from
                    for i in range(gene_pos_from,gene_pos_to+1):
                        genome_positions[gene_ann[(gene,i)]].add((gene,var))

                else:
                    #c.-16C>T
                    re_obj = re.match("c.(-[0-9]+)[A-Z]>[A-Z]", var)
                    gene_pos = int(re_obj.group(1))
                    genome_positions[gene_ann[(gene,gene_pos)]].add((gene,var))

            elif var[0]=="r":
                #rrl r.2814g>t
                re_obj = re.match("r.([0-9]+)[A-Za-z]+>[A-Za-z]+",var)
                gene_pos = int(re_obj.group(1))
                genome_positions[gene_ann[(gene,gene_pos)]].add((gene,var))


            elif var=="frameshift":
                pass
            elif var=="large_deletion":
                pass
            elif var[:19]=="any_missense_codon_":
                codon_pos = int(var.replace("any_missense_codon_",""))
                for pos in codon_ann[(gene,codon_pos)]:
                    genome_positions[pos].add((gene,var))
            else:
                print(gene,var)
                quit()

    return genome_positions


def rv2genes(bed_file):
    #Chromosome      759310  763325  Rv0667  rpoB    rifampicin
    rv2gene = {}
    for l in open(bed_file):
        row = l.strip().split()
        rv2gene[row[3]] = row[4]
    return rv2gene
