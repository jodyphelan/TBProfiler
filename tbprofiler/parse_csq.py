from __future__ import division
import sys
import re
import ann
import json

re_dp4 = re.compile("DP4=(\d+,\d+,\d+,\d+)")
re_mutation = re.compile("(\d+)(\w+)>\d+(\w+)")

def load_csq(self,pos_subset=None):
    #pos_subset = ["1","2","3"]
    obj_ann = ann.ann(self.params["annfile"],self.params["tabix"])
    results = {"loci":[],"variants":[]}
    for l in open(self.params["dr_vcffile"]):
        if l[0]=="#": continue
        arr = l.rstrip().split()
        if pos_subset and arr[1] not in pos_subset: continue
        dp4 = [int(x) for x in re_dp4.search(l).group(1).split(",")]
        dp4r = (dp4[2]+dp4[3])/(dp4[0]+dp4[1]+dp4[2]+dp4[3])

        var = {}
        var["pos"] = arr[1]
        var["chr"] = arr[0]
        var["ref"] = arr[3]
        var["alt"] = arr[4]
        var["dp4r"] = dp4r

        if "BCSQ" not in l:
            var["bcsq"] = "intergenic"
            results["variants"].append(var)
            results["loci"].append((arr[0],arr[1]))
            continue
        csq = re.search("BCSQ=(.*)",arr[7]).group(1).split("|")
        if csq[0][0]=="@": continue
        var["bcsq"] = csq[0]
        if csq[0]=="non_coding":
            var["csq"] = {"locus_tag":csq[1],"nt_change":"-","aa_change":"-"}
        elif csq[0]=="stop_lost&inframe_deletion":
            var["csq"] = {"locus_tag":csq[1],"nt_change":"-","aa_change":"-"}
        elif csq[0]=="inframe_deletion&start_lost":
            var["csq"] = {"locus_tag":csq[1],"nt_change":"-","aa_change":"-"}
        elif csq[0]=="frameshift":
            var["csq"] = {"locus_tag":csq[1],"nt_change":"-","aa_change":"-"}
        else:
            if len(csq)==7:
                var["csq"] = {"locus_tag":csq[1],"aa_change":csq[5],"nt_change":csq[6]}
            else:
                var["csq"] = {"locus_tag":csq[1],"nt_change":"-","aa_change":"-"}
        results["variants"].append(var)
        results["loci"].append((arr[0],arr[1]))

    if len(results["loci"])==0:
        print("Warning! No Variants")
    pos_tup = [("Chromosome",x["pos"]) for x in results["variants"]]
    dict_ann = obj_ann.pos2ann(self.params["stor_dir"],pos_tup)
    for var in results["variants"]:
        tann = dict_ann["Chromosome"][var["pos"]]
        if var["bcsq"]=="intergenic" or var["bcsq"]=="non_coding" or var["bcsq"]=="frameshift":
            locus_type = "intergenic" if tann["ncr"]=="inter" else "ncRNA"
            change = "%s%s>%s" % (tann["gene_nt"],var["ref"],var["alt"])
            var["csq"] = {"var_type":locus_type,"locus_tag":tann["rv"],"nt_change":change,"aa_change":"-"}
        var["csq"]["gene"] = tann["gene"]
    return results
