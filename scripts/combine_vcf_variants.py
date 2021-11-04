#! /usr/bin/env python3
from collections import defaultdict
import sys
import argparse
import re
import pathogenprofiler as pp

def gff_load_cds(gff):
    cds = []
    for l in open(gff):
        row = l.strip().split()
        if len(row)<3:
            continue
        if row[2]!="CDS":
            continue
        r = re.search("ID=(.+);+",l)
        id = r.group(1)
        cds.append({"chrom":row[0],"gene":id,"start":int(row[3]),"end":int(row[4]),"strand":row[6]})
    return cds

def get_overlapping_gene(chrom,pos,genes):
    genes = [g for g in genes if g["start"]<=pos and g["end"]>=pos and chrom==g["chrom"]]
    if len(genes)==0:
        return None
    else:
        return genes[0]

def get_codon_pos(chrom,pos,genes):
    g = get_overlapping_gene(chrom,pos,genes)
    if g==None:
        return (None,None)
    if g["strand"]=="+":
        codon_pos = (pos-g["start"])//3 + 1
    else:
        codon_pos = (g["end"] - pos )//3 + 1
    return (g["gene"],codon_pos)


def main(args):
    ref = pp.fasta(args.ref).fa_dict
    cds = gff_load_cds(args.gff)
    final_list = []
    coding = defaultdict(list)
    generator = pp.cmd_out(f"bcftools view {args.vcf}") if args.vcf else sys.stdin
    for l in generator:
        row = l.strip().split()
        if l[0]=="#":
            sys.stdout.write(l.strip()+"\n")
        elif len(row[3])>1 or len(row[4])>1:
            final_list.append(row)
        else:
            gene,cpos = get_codon_pos(row[0],int(row[1]),cds)
            if gene==None:
                final_list.append(row)
            else:
                coding[(gene,cpos)].append(row)
    

    for rows in coding.values():
        chrom = rows[0][0]
        pos = sorted([int(r[1]) for r in rows])

        ref_nucs = {p:ref[chrom][p-1] for p in range(pos[0],pos[-1]+1)}
        alt_nucs = ref_nucs.copy()
        for i,p in enumerate(pos):
            alt_nucs[p] = rows[i][4]
        new_row = rows[0]
        new_row[3] = "".join(ref_nucs.values())
        new_row[4] = "".join(alt_nucs.values())
        final_list.append(new_row)

    
    for row in sorted(final_list, key=lambda x:int(x[1])):
        sys.stdout.write("\t".join(row)+"\n")

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',type=str,help='')
parser.add_argument('--gff',type=str,help='',required = True)
parser.add_argument('--ref',type=str,help='',required = True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)