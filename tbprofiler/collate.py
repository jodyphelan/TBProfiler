import json
import files
from collections import defaultdict
import os
import re

_split = re.compile(r'[\0%s]' % re.escape(''.join(
    [os.path.sep, os.path.altsep or ''])))

def secure_filename(path):
    return _split.sub('', path)

class profiling_results:
    params = {}
    drugs = set()
    samples = []
    dr_drugs = {}
    def __init__(self,conf_file,samples_file,prefix,stor_dir,full_results,full_variant_results=False,db=None):
        tmp = json.load(open(conf_file))
        for x in tmp:
            self.params[x] = tmp[x]
        if db:
            self.params["dr_json"] = "%s.json" % db
            self.params["dr_bed_file"] = "%s.bed" % db
            files.filecheck(self.params["dr_json"])
            files.filecheck(self.params["dr_bed_file"])
        for l in open(self.params["dr_bed_file"]):
            arr = l.rstrip().split()
            for d in arr[5].split(";"):
                self.drugs.add(d)
        for s in open(samples_file):
            self.samples.append(s.rstrip())
        self.prefix = prefix
        self.stor_dir = stor_dir
        self.full_results = full_results
        self.full_variant_results = full_variant_results

    def results2tab(self):
        results = defaultdict(dict)
        linresults = defaultdict(dict)
        dr_variants = defaultdict(lambda:defaultdict(dict))
        dr_variants_set = set()

        for s in self.samples:
            for d in self.drugs:
                results[s][d] = set()
        for s in self.samples:
            temp = json.load(open("%s/results/%s.results.json" % (self.stor_dir,s)))
            # p = "%s/results\0/%s.results.json" % (self.stor_dir,s)
            # p = re.sub(r'\0', '', p)
            # temp = json.load(open(p))
            for x in temp["small_variants_dr"]:
                for d in x["drug"].split(";"):
                    dr_variants[x["gene"]][x["change"]][s] = x["freq"]
                    dr_variants_set.add((x["gene"],x["change"]))
                    results[s][d].add("%s_%s" % (x["gene"],x["change"]) if self.full_results else "R")
            for x in temp["del"]:
                for d in x["drug"].split(";"):
                    results[s][d].add("large_deletion_%s" % x["gene"] if self.full_results else "R")
            for d in self.drugs:
                results[s][d] = ", ".join(results[s][d]) if len(results[s][d])>0 else "-"
            def get_dots(x):
                return len([c for c in x if c=="."])
            def derive_path(x):
                return [".".join(x.split(".")[:i])for i in range(1,len(x.split(".")))] + [x]
            lin = temp["lineage"]
            lin_fracs = {}
            pool = []
            for l in lin:
                pool.append(l["lin"])
                lin_fracs[l["lin"]] = float(l["frac"])

            tree = {}
            while len(pool)>0:
                pool = sorted(pool,key=lambda x:get_dots(x))
                node = pool.pop(0)
                path = derive_path(node)
                parent = path[-2] if len(path)>1 else None
                if not parent:
                    tree[node] = {}
                else:
                    tmp = tree
                    for i in range(len(path)-1):
                        if path[i] not in tmp:
                            break
                        tmp = tmp[path[i]]
                    tmp[node] = {}
            linresults[s]["main"] = sorted(tree.keys(),key=lambda x:lin_fracs[x],reverse=True)[0]
            tmp = tree[linresults[s]["main"]]
            n = "NA"
            while len(tmp.keys())>0:
                n = sorted(tmp.keys(),key=lambda x:lin_fracs[x])[0]
                tmp = tmp[n]
            linresults[s]["sublin"] = n
            dr_drugs = [x["drug"] for x in temp["small_variants_dr"]]
            self.dr_drugs[s] = dr_drugs
            MDR = "R" if ("ISONIAZID" in dr_drugs and "RIFAMPICIN" in dr_drugs) else "-"
            XDR = "R" if MDR=="R" and ( "AMIKACIN" in dr_drugs or "KANAMYCIN" in dr_drugs or "CAPREOMYCIN" in dr_drugs ) and ( "FLUOROQUINOLONES" in dr_drugs) else "-"
            drtype = "Sensitive"
            if XDR=="R":
                drtype="XDR"
            elif MDR=="R":
                drtype="MDR"
            elif len(dr_drugs)>0:
                drtype="Drug-resistant"
            results[s]["XDR"] = XDR
            results[s]["MDR"] = MDR
            results[s]["drtype"] = drtype

        if self.full_variant_results:
            all_vars = json.load(open(self.params["dr_json"]))
            lt2gene = {}
            for l in open(self.params["dr_bed_file"]):
                #Chromosome      5240    7267    Rv0005  gyrB    FLUOROQUINOLONES
                row = l.rstrip().split()
                lt2gene[row[3]] = row[4] if row[4]!="." else row[3]
            dr_variants_set = set()
            for g in all_vars:
                for c in all_vars[g]:
                    dr_variants_set.add((lt2gene[g],c))
        VAR = open(self.prefix+".variants.txt","w")
        VAR.write("sample\t%s\n" % ("\t".join(["%s_%s" % (g,c) for g,c in sorted(dr_variants_set,key=lambda x: x[0])])))
        for s in self.samples:
            VAR.write("%s\t%s\n" % (s,"\t".join(["%.3f" % dr_variants[gene][change][s] if gene in dr_variants and change in dr_variants[gene] and s in dr_variants[gene][change] else "0" for gene,change in sorted(dr_variants_set,key=lambda x: x[0])])))
        VAR.close()
        o = open(self.prefix+".txt","w")
        o.write("sample\tmain_lineage\tsub_lineage\tDR_type\tMDR\tXDR\t%s" % "\t".join(self.drugs)+"\n")
        for s in self.samples:
            o.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(s,linresults[s]["main"],linresults[s]["sublin"],results[s]["drtype"],results[s]["MDR"],results[s]["XDR"],"\t".join([results[s][x] for x in self.drugs])))
        o.close()
        json.dump(results,open(self.prefix+".json","w"))
        lineage_cols = {"lineage1":"#104577","lineage2":"#ab2323","lineage3":"#18a68c","lineage4":"#f68e51","lineage5":"#7cb5d2","lineage6":"#fde05e","lineage7":"#bc94b7","lineageBOV":"#f8e0c8","lineageOther":"#000000"}
        o = open(self.prefix+".lineage.itol.txt","w")
        o.write("""DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL    Lineage
COLOR    #ff0000

LEGEND_TITLE    Lineage
LEGEND_SHAPES    1    1    1    1    1    1    1    1    1
LEGEND_COLORS    %(lineage1)s    %(lineage2)s    %(lineage3)s    %(lineage4)s    %(lineage5)s    %(lineage6)s    %(lineage7)s    %(lineageBOV)s    %(lineageOther)s
LEGEND_LABELS    Lineage1    Lineage2    Lineage3    Lineage4    Lineage5    Lineage6    Lineage7    Bovis    Other

DATA
""" % lineage_cols)
        for s in self.samples:
            o.write("%s\t%s\n" % (s,lineage_cols.get(linresults[s]["main"],"#000000")))
        o.close()

        o = open(self.prefix+".dr.itol.txt","w")
        dr_cols = {"Sensitive":"#80FF00","Drug-resistant":"#00FFFF","MDR":"#8000FF","XDR":"#FF0000"}
        o.write("""DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL    Drug-Resistance
COLOR    #ff0000

LEGEND_TITLE    Drug resistance
LEGEND_SHAPES    1    1    1    1
LEGEND_COLORS    #80FF00    #00FFFF    #8000FF    #FF0000
LEGEND_LABELS    Sensitive    Drug-resisant    MDR    XDR

DATA
""")
        for s in self.samples:
            o.write("%s\t%s\n" % (s,dr_cols.get(results[s]["drtype"],"#000000")))
        o.close()

        self.drugs = ["RIFAMPICIN","ISONIAZID","ETHAMBUTOL","PYRAZINAMIDE","STREPTOMYCIN","FLUOROQUINOLONES","AMINOGLYCOSIDES","KANAMYCIN","AMIKACIN","CAPREOMYCIN","ETHIONAMIDE","PARA-AMINOSALISYLIC_ACID","CLOFAZIMINE","LINEZOLID","BEDAQUILINE"]
        o = open(self.prefix+".dr.indiv.itol.txt","w")
        dr_cols = {"Sensitive":"#80FF00","Drug-resistant":"#00FFFF","MDR":"#8000FF","XDR":"#FF0000"}
        legend_shapes = "\t".join(["2" for x in self.drugs])
        legend_colours = "\t".join(["black" for x in self.drugs])
        legend_labels = "\t".join(self.drugs)
        o.write("""DATASET_BINARY
SEPARATOR TAB
DATASET_LABEL    Drugs
COLOR    #ff0000

SHOW_LABELS    1
FIELD_SHAPES    %s
FIELD_COLORS    %s
FIELD_LABELS    %s

DATA
""" % (legend_shapes,legend_colours,legend_labels))
        for s in self.samples:
            o.write("%s\t%s\n" % (s,"\t".join(["1" if d in self.dr_drugs[s] else "0" for d in self.drugs])))
        o.close()
