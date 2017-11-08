import sys
import json
import re
from collections import defaultdict
sys.path.insert(0,"../tbprofiler/")
from tbprofiler import ann

amino_acids = ['Cys', 'Ile', 'Ser', 'Val', 'Gly', 'Gln', 'Pro', 'Lys', 'Stop', 'Thr', 'Phe', 'Ala', 'Met', 'Asp', 'His', 'Leu', 'Arg', 'Trp', 'Glu', 'Asn', 'Tyr']
aa_long2short = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V","Stop":"*"}
re_mut = re.compile("([A-Za-z]+)([\s-]*[0-9]+)([A-Za-z]+)")
infile = sys.argv[1]
prefix = sys.argv[2]
json_db = "%s.json" % prefix
bedfile = "%s.bed" % prefix

drdb = defaultdict(lambda:defaultdict(list))

positions = set()
for l in open(infile):
	arr = l.rstrip().split()
	for p in arr[1].split("/"):
		if p=="-": continue
		positions.add(int(p))

obj_ann = ann.ann("/Users/jody/github/TBProfiler/ref/MTB-h37rv_asm19595v2-eg18.tab.ann.gz","/Users/jody/github/TBProfiler/bin/tabix")

annot = obj_ann.pos2ann([("Chromosome",str(x)) for x in sorted(positions)])
loci = {}
for l in open(infile):
	#AMINOGLYCOSIDES 1473247 C	   A	   rrs	 C1402A
	arr = l.rstrip().split()
	re_obj = re_mut.match(arr[5])
	ref = re_obj.group(1)
	gene_pos = re_obj.group(2)
	alt = re_obj.group(3)
	if arr[1]=="-": continue
	tmp = annot["Chromosome"][arr[1].split("/")[0]]
	drugs = loci[tmp["rv"]][4] if tmp["rv"] in loci else set()
	drugs.add(arr[0])
	loci[tmp["rv"]] = (tmp["start"],tmp["end"],tmp["rv"],tmp["gene"],drugs)
	if ref in amino_acids:
		mutation = "%s%s>%s%s" % (gene_pos,aa_long2short[ref],gene_pos,aa_long2short[alt])
		drdb[tmp["rv"]][mutation].append(arr[0])
	else:
		mutation = "%s%s>%s" % (tmp["gene_nt"],arr[2],arr[3])
		drdb[tmp["rv"]][mutation].append(arr[0])
json.dump(drdb,open(json_db,"w"))
O = open(bedfile,"w")
for rv in sorted(loci,key=lambda x:int(loci[x][0])):
	tmp = loci[rv]
	O.write("Chromosome\t%s\t%s\t%s\t%s\t%s\n" % (tmp[0],tmp[1],tmp[2],tmp[3],";".join(tmp[4])))
