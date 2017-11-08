import json


def db_compare(self,res):
    drugs = set()
    variants = {"dr":[],"other":[]}
    drdb = json.load(open(self.params["dr_json"]))
    for var in res["variants"]:
        tmp = "aa_change" if var["csq"]["aa_change"] != "-" else "nt_change"
        if var["csq"]["locus_tag"] in drdb:
            if var["csq"]["aa_change"] in drdb[var["csq"]["locus_tag"]] or var["csq"]["nt_change"] in drdb[var["csq"]["locus_tag"]]:
                #tmp = "aa_change" if var["csq"]["aa_change"] in drdb[var["csq"]["locus_tag"]] else "nt_change"
                var["drugs"] = drdb[var["csq"]["locus_tag"]][var["csq"][tmp]]
                #variants["dr"].append(var)
                for d in drdb[var["csq"]["locus_tag"]][var["csq"][tmp]]:
                    variants["dr"].append({"locus_tag":var["csq"]["locus_tag"],"chr":"Chromosome","change":var["csq"][tmp],"freq":var["dp4r"],"gene":var["csq"]["gene"],"type":var["bcsq"],"genome_pos":var["pos"],"drug":d})
                    drugs.add(d)
            else:
                variants["other"].append({"locus_tag":var["csq"]["locus_tag"],"chr":"Chromosome","change":var["csq"][tmp],"freq":var["dp4r"],"gene":var["csq"]["gene"],"type":var["bcsq"],"genome_pos":var["pos"]})
    return variants
