import sys
import json
import re
from collections import defaultdict
import ann

amino_acids = ['Cys', 'Ile', 'Ser', 'Val', 'Gly', 'Gln', 'Pro', 'Lys', 'Stop', 'Thr', 'Phe', 'Ala', 'Met', 'Asp', 'His', 'Leu', 'Arg', 'Trp', 'Glu', 'Asn', 'Tyr']
aa_long2short = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V","Stop":"*"}
re_mut = re.compile("([A-Za-z]+)([\s-]*[0-9]+)([A-Za-z]+)")
infile = sys.argv[1]
outfile = sys.argv[2]
drdb = defaultdict(lambda:defaultdict(list))

positions = set()
for l in open(infile):
    arr = l.rstrip().split()
    for p in arr[1].split("/"):
        if p=="-": continue
        positions.add(int(p))

obj_ann = ann.ann("/Users/jody/github/TBProfiler/ref/MTB-h37rv_asm19595v2-eg18.tab.ann.gz","/Users/jody/github/TBProfiler/bin/tabix")

annot = obj_ann.pos2ann([("Chromosome",str(x)) for x in sorted(positions)])

for l in open(infile):
#    print l.rstrip()
    #AMINOGLYCOSIDES 1473247 C       A       rrs     C1402A
    arr = l.rstrip().split()
    re_obj = re_mut.match(arr[5])
    ref = re_obj.group(1)
    gene_pos = re_obj.group(2)
    alt = re_obj.group(3)
    if arr[1]=="-": continue
#    print "%s\t%s\t%s" % (ref,gene_pos,alt)
    tmp = annot["Chromosome"][arr[1].split("/")[0]]
    if ref in amino_acids:
        mutation = "%s%s>%s%s" % (gene_pos,aa_long2short[ref],gene_pos,aa_long2short[alt])
        drdb[tmp["rv"]][mutation].append(arr[0])
    else:
        mutation = "%s%s>%s" % (tmp["gene_nt"],arr[2],arr[3])
        drdb[tmp["rv"]][mutation].append(arr[0])
json.dump(drdb,open(outfile,"w"))
