from __future__ import division
from files import *
from collections import defaultdict
def deletions(self,min_gene_frac=0.9):

    lt_bed = defaultdict(dict)
    lt_cov = defaultdict(lambda : defaultdict(int))
    lt_drugs = {}
    lt2gene = {}
    for l in open(self.params["dr_bed_file"]):
        chrom,start,end,lt,gene,drugs = l.rstrip().split()
        lt_drugs[lt] = drugs
        for pos in range(int(start),int(end)+1):
            lt_cov[lt][str(pos)] = 0
            lt_bed[chrom][str(pos)] = lt
            lt2gene[lt] = gene

    cmd = "%(samtools)s depth -b %(dr_bed_file)s %(bamfile)s > %(depthfile)s" % self.params
    run_cmd(cmd,verbose=self.params["verbose"])

    for l in open(self.params["depthfile"]):
        chrom,pos,depth =  l.rstrip().split()
        lt_cov[lt_bed[chrom][pos]][pos] = depth
    results = []
    for lt in lt_cov:
        arr_cov = lt_cov[lt].values()
        int_good_cov = len(filter(lambda x: x>3,arr_cov))
        float_gene_present = int_good_cov/len(arr_cov)

        if float_gene_present<min_gene_frac:
            if lt=="Rv0667":
                print "Missing coverage on essential gene rpoB (%s). Stopping pipeline. Check for contamination" % float_gene_present
                quit()
                return []
            results.append({"drug":lt_drugs[lt],"gene":lt2gene[lt],"locus_tag":lt,"len":len(arr_cov),"good_cov":int_good_cov,"float_gene_present":float_gene_present})

    return results
