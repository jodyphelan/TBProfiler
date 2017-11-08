import sys
import re
i = 0
for l in open(sys.argv[1]):
    #
    if l[0]=="#": print l.strip(); continue
    arr = l.strip().split()
    if arr[2] != "gene": continue
    gene_no = "gene%s" % i
    i+=1
    rv =  re.search("locus_tag=(\w+)",l).group(1)
    gene_name = re.search("Name=(\w+)",l).group(1) if re.search("Name=(\w+)",l)!=None else rv
    gene_biotype = re.search("gene_biotype=(\w+)",l).group(1)
    print "Chromosome\t.\tgene\t%s\t%s\t.\t%s\t.\tID=gene:%s;Name=%s;gene=%s;locus_tag=%s;gene_biotype=%s;" % (arr[3],arr[4],arr[6],gene_no,rv,gene_name,rv,gene_biotype)
    print "Chromosome\t.\ttranscript\t%s\t%s\t.\t%s\t.\tID=transcript:%s;Parent=gene:%s;gene=%s;locus_tag=%s;gene_biotype=%s;" % (arr[3],arr[4],arr[6],gene_no,gene_no,gene_name,rv,gene_biotype)
    print "Chromosome\t.\tCDS\t%s\t%s\t.\t%s\t.\tParent=transcript:%s;gene=%s" % (arr[3],arr[4],arr[6],gene_no,gene_name)
#Chromosome	RefSeq	gene	7302	9818	.	+	.	ID=gene:gene5;Dbxref=GeneID:887105;Name=gyrA;experiment=DESCRIPTION:Mutation analysis[PMID:8031045];gbkey=Gene;gene=gyrA;gene_biotype=protein_coding;locus_tag=Rv0006
#Chromosome	RefSeq	transcript	7302	9818	.	+	.	ID=transcript:tr5;Parent=gene:gene5;Dbxref=GeneID:887105;Name=gyrA;experiment=DESCRIPTION:Mutation analysis[PMID:8031045];gbkey=Gene;gene=gyrA;gene_biotype=protein_coding;locus_tag=Rv0006
#Chromosome	RefSeq	CDS	7302	9818	.	+	0	ID=cds5;Parent=transcript:tr5;Dbxref=Genbank:NP_214520.1,GeneID:887105;Name=NP_214520.1;experiment=COORDINATES:Mass spectrometry[PMID:21969609];gbkey=CDS;gene=gyrA;inference=protein motif:PROSITE:PS00018;product=DNA gyrase subunit A;protein_id=NP_214520.1;transl_table=11
