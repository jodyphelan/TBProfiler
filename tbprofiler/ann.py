import sys
import random
import subprocess
import os
from collections import defaultdict
r = random.SystemRandom()

class ann:
    annfile = ""
    tabix = ""
    def __init__(self,filename,tabix):
        self.annfile = filename
        self.tabix = tabix
    def pos2ann(self,stor_dir, pos_tuple_list):
        #pos_tuple_list = [("Chromosome",1),("Chromosome",2)]
        if len(pos_tuple_list)<5000:
            p1offset = -1
            p2offset = 0
            stype = "R"
        else:
            p1offset = 0
            p2offset = 1
            stype = "T"

        num = r.randint(1,1000000)+r.randint(1,1000000)
        temp_bed_file = "%s/temp.%s.bed" % (stor_dir,num)
        OUT = open(temp_bed_file,"w")
        for chrom,pos in pos_tuple_list:
            OUT.write("%s\t%s\t%s\n" % (chrom,int(pos)+p1offset,int(pos)+p2offset))
        OUT.close()
        results = defaultdict(dict)
        for l in subprocess.Popen("%s %s -%s %s " % (self.tabix,self.annfile,stype,temp_bed_file),stdout=subprocess.PIPE,shell=True).stdout:
            # added .decode() to convert binary output of subprocess.Popen to string
            l = l.decode()
            arr = l.rstrip().split()
            results[arr[0]][arr[1]] = {"change_pos":arr[7],"ref_nt":arr[2],"ref_codon":arr[6],"ref_aa":arr[11],"chr":arr[0],"pos":arr[1],"rv":arr[15],"gene":arr[16],"gene_syn":arr[17],"ncr":arr[18],"start":arr[19],"end":arr[20],"strand":arr[21],"codon_num":arr[24],"gene_nt":arr[25],"operon":arr[26]}
            
        os.remove(temp_bed_file)
        return results
