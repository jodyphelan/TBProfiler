#! /usr/bin/python
import sys
import re
import subprocess
from collections import defaultdict
import json

tabix = "/Users/jody/github/TBProfiler/bin/tabix"
ref_annotation = "/Users/jody/github/TBProfiler/ref/MTB-h37rv_asm19595v2-eg18.tab.ann.gz"
chrom = "Chromosome"
aa_convert = {"Ile":"I","Leu":"L","Val":"V","Phe":"F","Met":"M","Cys":"C","Ala":"A","Gly":"G","Pro":"P","Thr":"T","Ser":"S","Tyr":"Y","Trp":"W","Gln":"Q","Asn":"N","His":"H","Glu":"E","Asp":"D","Lys":"K","Arg":"R","Stop":"*"}
rc = {"A":"T","C":"G","G":"C","T":"A"}
amino_acids = ['Cys', 'Ile', 'Ser', 'Val', 'Gly', 'Gln', 'Pro', 'Lys', 'Stop', 'Thr', 'Phe', 'Ala', 'Met', 'Asp', 'His', 'Leu', 'Arg', 'Trp', 'Glu', 'Asn', 'Tyr']
aa2codons = {'Cys': ['TGT', 'TGC'], 'Asp': ['GAT', 'GAC'], 'Ile': ['ATC', 'ATA', 'ATT'], 'Ser': ['AGC', 'AGT', 'TCT', 'TCG', 'TCC', 'TCA'], 'Gln': ['CAA', 'CAG'], 'Lys': ['AAG', 'AAA'], 'Pro': ['CCT', 'CCG', 'CCA', 'CCC'], 'Stop': ['TAG', 'TAA', 'TGA'], 'Thr': ['ACA', 'ACG', 'ACT', 'ACC'], 'Phe': ['TTT', 'TTC'], 'Ala': ['GCA', 'GCC', 'GCG', 'GCT'], 'Met': ['ATG'], 'Gly': ['GGT', 'GGG', 'GGA', 'GGC'], 'His': ['CAT', 'CAC'], 'Leu': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'], 'Arg': ['AGG', 'AGA', 'CGA', 'CGG', 'CGT', 'CGC'], 'Trp': ['TGG'], 'Val': ['GTA', 'GTC', 'GTG', 'GTT'], 'Asn': ['AAC', 'AAT'], 'Tyr': ['TAT', 'TAC'], 'Glu': ['GAG', 'GAA']}

def recode_indel(ref,alt):
	if len(ref)<len(alt):
		ldiff = len(alt)-len(ref)
		return "%s+%s%s" % (ref[:1],ldiff,alt[1:1+ldiff])
	else:
		ldiff = len(ref)-len(alt)
		return "%s-%s%s" % (ref[:1],ldiff,ref[1:1+ldiff])


def load_ann(bed_name):
	ann_dict = defaultdict(dict)
	annCMD = "%s %s -T %s" % (tabix,ref_annotation,bed_name)
	annPIPE = subprocess.Popen(annCMD,shell=True,stdout=subprocess.PIPE)
	for l in annPIPE.stdout:
		arr = l.rstrip().split()
		nuc_obj = {}
		nuc_obj[arr[2]] = {"codon":arr[6],"aa":arr[11]}
		nuc_obj[arr[3]] = {"codon":arr[8],"aa":arr[12]}
		nuc_obj[arr[4]] = {"codon":arr[9],"aa":arr[13]}
		nuc_obj[arr[5]] = {"codon":arr[10],"aa":arr[14]}
		ann_dict[arr[0]][arr[1]] = {"ref_nt":arr[2],"ref_codon":arr[6],"ref_aa":arr[11],"chr":arr[0],"pos":arr[1],"ann":nuc_obj,"rv":arr[15],"gene":arr[16],"gene_syn":arr[17],"ncr":arr[18],"start":arr[19],"end":arr[20],"strand":arr[21],"drug":arr[22],"ppe":arr[23],"codon_num":arr[24],"gene_nt":arr[25],"operon":arr[26]}
	return ann_dict



bed_pos = set()

def revcom(x):
	return "".join([{"A":"T","C":"G","G":"C","T":"A"}[x[i]] for i in range(len(x)-1,-1,-1)])

ann = load_ann("temp.bed")

temp_drdb = defaultdict(list)

for l in open("drdb.txt"):
	print l
	#ETHIONAMIDE	1674481	T	G	inhA	Ser94Ala
	#ETHIONAMIDE	1674484/1674485	AT	CC	inhA	Ile95Pro
	arr = l.rstrip().split()
	if arr[2]=="-": continue
	positions = arr[1].split("/")
	for p in positions:
		bed_pos.add((ann[chrom][p]["start"],ann[chrom][p]["end"],ann[chrom][p]["rv"],ann[chrom][p]["gene"]))
	pos = positions[0]
	t_ann = ann["Chromosome"][pos]
	ncr = ann["Chromosome"][pos]["ncr"]
	ref_codon = ann[chrom][pos]["ref_codon"]
	strand = ann[chrom][pos]["strand"]
	vartype = "SNP" if len(arr[2])==1 and len(arr[3])==1 else "INDEL"
	gene = ann[chrom][pos]["gene"]
	rv = ann[chrom][pos]["rv"]
	unit_pos =  []


	if ann[chrom][pos]["codon_num"]!=".":
		unit_pos = [ann[chrom][str(ipos)]["pos"] for ipos in range(int(ann[chrom][pos]["start"]),int(ann[chrom][pos]["end"])+1) if ann[chrom][str(ipos)]["codon_num"]==ann[chrom][pos]["codon_num"]]
	else:
		unit_pos = [ann[chrom][pos]["gene_nt"]]

	re_obj = re.search("([A-Za-z]+)(-*\d+)([A-Za-z]+)",l)
	ref,change_pos,alt = re_obj.group(1),re_obj.group(2),re_obj.group(3)

	print ref_codon
	if ann[chrom][pos]["strand"]=="-":
		ref_codon = revcom(ref_codon)
	print ref_codon
	key = ""

	if ref in amino_acids and alt in amino_acids:
		if ref==alt:
			change_str = "%s%s%s" % (ref,change_pos,alt)
			key = json.dumps({"genome_pos":[pos],"ref_nt":[ann[chrom][pos]["ref_nt"]],"alt_nt":[arr[3]],"locus_tag":rv,"gene":gene,"chr":"Chromosome","change":change_str})
		else:
			for alt_codon in aa2codons[alt]:
				if ann[chrom][pos]["strand"]=="-":
					alt_codon = revcom(alt_codon)
				print ref_codon
				alt_pos = [unit_pos[j] for j in [i for i in range(3) if ref_codon[i]!=alt_codon[i]]]
				ref_nts = [ref_codon[i] for i in range(3) if ref_codon[i]!=alt_codon[i]]
				alt_nts = [alt_codon[i] for i in range(3) if ref_codon[i]!=alt_codon[i]]
				print ref_nts
				change_str = "%s%s%s" % (ref,change_pos,alt)
				key = json.dumps({"genome_pos":alt_pos,"ref_nt":ref_nts,"alt_nt":alt_nts,"locus_tag":rv,"gene":gene,"chr":"Chromosome","change":change_str})
	elif len(alt)!=len(ref):
		change_str = "%s%s%s" % (ann[chrom][pos]["ref_nt"],ann[chrom][pos]["gene_nt"],recode_indel(ref,alt))
		key = json.dumps({"genome_pos":[pos],"ref_nt":[ann[chrom][pos]["ref_nt"]],"alt_nt":[recode_indel(ref,alt)],"locus_tag":rv,"gene":gene,"chr":"Chromosome","change":change_str})

	else:
		change_str = "%s%s%s" % (ref,change_pos,alt)
 		key = json.dumps({"genome_pos":[pos],"ref_nt":[ann[chrom][pos]["ref_nt"]],"alt_nt":[arr[3]],"locus_tag":rv,"gene":gene,"chr":"Chromosome","change":change_str})
	temp_drdb[key].append(arr[0])

drdb = []
for key in sorted(temp_drdb,key=lambda x: int(json.loads(x)["genome_pos"][0])):
	obj = json.loads(key)
	obj["drug"] = temp_drdb[key]
	drdb.append(obj)
json.dump(drdb,open("drdb.json","w"))
with open("drdb.bed","w") as o:
	for p in sorted(list(bed_pos),key=lambda x:int(x[0])):
		o.write("Chromosome\t%s\t%s\t%s\t%s\n" % (p[0],p[1],p[2],p[3]))
