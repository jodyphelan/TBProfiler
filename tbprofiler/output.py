import json

def load_drgus(dr_bed_file):
    drugs = set()
    for l in open(dr_bed_file):
        arr = l.rstrip().split()
        for d in arr[5].split(";"):
            drugs.add(d)
    return sorted(list(drugs))


def write_results_2(self):
    drugs = load_drgus(self.params["dr_bed_file"])
    o = open(self.params["txt_results"],"w")
    o.write("#### %(prefix)s results ####\n\n" % self.params)
    o.write("#### Drug resistance associated small variants ####\n")
    for d in drugs:
        dr_vars = [x for x in self.small_dr_variants["dr"] if x["drug"]==d]
        res = "R" if len(dr_vars)>0 else ""
        muts = ["%(gene)s (%(change)s, %(freq)s)" % x for x in dr_vars]
        o.write("%s\t%s\t%s\n" % (d,res,", ".join(muts)))
    dr_drugs = [x["drug"] for x in self.small_dr_variants["dr"]]

    MDR = "R" if ("ISONIAZID" in dr_drugs and "RIFAMPICIN" in dr_drugs) else ""
    XDR = "R" if MDR=="R" and ( "AMIKACIN" in dr_drugs or "KANAMYCIN" in dr_drugs or "CAPREOMYCIN" in dr_drugs ) and ( "FLUOROQUINOLONES" in dr_drugs) else ""
    o.write("MDR\t%s\n" % MDR)
    o.write("XDR\t%s\n" % XDR)

    o.write("\n\n#### Big Deletions in candidate genes ####\n")
    for l in self.deletions:
        o.write("\t".join([str(l[x]) for x in ["drug","gene","len","good_cov",'float_gene_present']])+"\n")

    o.write("\n\n#### Lineage ####\n")
    o.write("\n".join(["%s\t%s\t%s\t%s\t%s" % (x["lin"],x["frac"],x["family"],x["spoligotype"],x["rd"]) for x in self.lineage]))

    o.write("\n\n#### Other small variants in candidate genes ####\n")
    for l in self.small_dr_variants["other"]:
        o.write("\t".join([str(l[x]) for x in ["chr","genome_pos","gene","locus_tag","change","type","freq"]])+"\n")
    o.close()

    json_dict = {"small_variants_other":self.small_dr_variants["other"],"small_variants_dr":self.small_dr_variants["dr"],"del":self.deletions,"lineage":self.lineage,"id":self.params["prefix"]}

    json.dump(json_dict,open(self.params["json_results"],"w"))

def write_results_1(self):
    drugs = load_drgus(self.params["dr_bed_file"])
    o = open(self.params["txt_results"],"w")
    o.write("#### %(prefix)s results ####\n\n" % self.params)
    o.write("#### Drug resistance associated small variants ####\n")
    for l in self.small_dr_variants["dr"]:
        o.write("\t".join([str(l[x]) for x in ["drug","chr","genome_pos","gene","locus_tag","change","type","freq"]])+"\n")
    o.write("\n\n#### Big Deletions in candidate genes ####\n")
    for l in self.deletions:
        o.write("\t".join([str(l[x]) for x in ["drug","gene","len","good_cov",'float_gene_present']])+"\n")
    o.write("\n\n#### Lineage ####\n")
    o.write("\n".join(["%s\t%s\t%s\t%s\t%s" % (x["lin"],x["frac"],x["family"],x["spoligotype"],x["rd"]) for x in self.lineage]))
    o.write("\n\n#### Other small variants in candidate genes ####\n")
    for l in self.small_dr_variants["other"]:
        o.write("\t".join([str(l[x]) for x in ["chr","genome_pos","gene","locus_tag","change","type","freq"]])+"\n")
    o.close()

    json_dict = {"small_variants_other":self.small_dr_variants["other"],"small_variants_dr":self.small_dr_variants["dr"],"del":self.deletions,"lineage":self.lineage,"id":self.params["prefix"]}

    json.dump(json_dict,open(self.params["json_results"],"w"))
