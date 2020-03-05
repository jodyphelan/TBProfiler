import json
import os
import sys
from collections import defaultdict
from tqdm import tqdm

def collate_results(prefix,conf,dir="./results",sample_file=None,full_results=True,full_variant_results=True,reporting_af=0.1):
    if not os.path.isdir(dir):
        sys.stderr.write("\nERROR: Can't find directory %s\n" % dir )
        exit()
    set_all_drugs = set()
    for l in open(conf["bed"]):
        arr = l.rstrip().split()
        for d in arr[5].split(","):
            set_all_drugs.add(d)
    tmp_drugs = ["rifampicin", "isoniazid", "pyrazinamide", "ethambutol", "streptomycin", "fluoroquinolones", "moxifloxacin", "ofloxacin", "levofloxacin", "ciprofloxacin", "aminoglycosides", "amikacin", "kanamycin", "capreomycin", "ethionamide", "para-aminosalicylic_acid", "cycloserine", "linezolid"]
    drug_list = []
    for d in tmp_drugs:
        if d in set_all_drugs:
            drug_list.append(d)
    for d in sorted(list(set_all_drugs)):
        if d not in drug_list:
            drug_list.append(d)

    if sample_file:
        samples = [x.rstrip() for x in open(sample_file).readlines()]
    else:
        samples = [x.replace(".results.json","") for x in os.listdir("%s/" % dir) if x[-13:]==".results.json"]

    results = defaultdict(dict)
    dr_variants = defaultdict(lambda:defaultdict(dict))
    dr_variants_set = set()
    dr_drugs = {}
    for s in samples:
        for d in drug_list:
            results[s][d] = set()
    for s in tqdm(samples):
        temp = json.load(open("%s/%s.results.json" % (dir,s)))
        for x in temp["dr_variants"]:
            if x["freq"]>reporting_af:
                dr_variants[x["gene"]][x["change"]][s] = x["freq"]
                dr_variants_set.add((x["gene"],x["change"]))
                results[s][x["drug"]].add("%s_%s" % (x["gene"],x["change"]) if full_results else "R")
        for d in drug_list:
            results[s][d] = ", ".join(results[s][d]) if len(results[s][d])>0 else "-"
            results[s]["main_lin"] = temp["main_lin"]
            results[s]["sublin"] = temp["sublin"]
            results[s]["drtype"] = temp["drtype"]
            results[s]["MDR"] = temp["MDR"]
            results[s]["XDR"] = temp["XDR"]
        dr_drugs[s] = [x["drug"] for x in temp["dr_variants"]]
    if full_variant_results:

        all_vars = json.load(open(conf["json_db"]))
        lt2gene = {}
        for l in open(conf["bed"]):
            #Chromosome      5240    7267    Rv0005  gyrB    FLUOROQUINOLONES
            row = l.rstrip().split()
            lt2gene[row[3]] = row[4] if row[4]!="." else row[3]
        for gene in all_vars:
            for mutation in all_vars[gene]:
                dr_variants_set.add((lt2gene[gene],all_vars[gene][mutation]["hgvs_mutation"]))
        VAR = open(prefix+".variants.txt","w")
        VAR.write("sample\t%s\n" % ("\t".join(["%s_%s" % (g,c) for g,c in sorted(dr_variants_set,key=lambda x: x[0])])))
        for s in samples:
            VAR.write("%s\t%s\n" % (s,"\t".join(["%.3f" % dr_variants[gene][change][s] if gene in dr_variants and change in dr_variants[gene] and s in dr_variants[gene][change] else "0" for gene,change in sorted(dr_variants_set,key=lambda x: x[0])])))
        VAR.close()

    OUT = open(prefix+".txt","w")
    OUT.write("sample\tmain_lineage\tsub_lineage\tDR_type\tMDR\tXDR\t%s" % "\t".join(drug_list)+"\n")
    for s in samples:
        OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(s,results[s]["main_lin"],results[s]["sublin"],results[s]["drtype"],results[s]["MDR"],results[s]["XDR"],"\t".join([results[s][x] for x in drug_list])))
    OUT.close()
    json.dump(results,open(prefix+".json","w"))
    lineage_cols = {"lineage1":"#104577","lineage2":"#ab2323","lineage3":"#18a68c","lineage4":"#f68e51","lineage5":"#7cb5d2","lineage6":"#fde05e","lineage7":"#bc94b7","lineageBOV":"#f8e0c8","lineageOther":"#000000"}
    OUT = open(prefix+".lineage.itol.txt","w")
    OUT.write("""DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL\tLineage
COLOR\t#ff0000

LEGEND_TITLE\tLineage
LEGEND_SHAPES\t1\t1\t1\t1\t1\t1\t1\t1\t1
LEGEND_COLORS\t%(lineage1)s\t%(lineage2)s\t%(lineage3)s\t%(lineage4)s\t%(lineage5)s\t%(lineage6)s\t%(lineage7)s\t%(lineageBOV)s\t%(lineageOther)s
LEGEND_LABELS\tLineage1\tLineage2\tLineage3\tLineage4\tLineage5\tLineage6\tLineage7\tBovis\tOther

DATA
""" % lineage_cols)
    for s in samples:
        OUT.write("%s\t%s\n" % (s,lineage_cols.get(results[s]["main_lin"],"#000000")))
    OUT.close()

    OUT = open(prefix+".dr.itol.txt","w")
    dr_cols = {"Sensitive":"#80FF00","Drug-resistant":"#00FFFF","MDR":"#8000FF","XDR":"#FF0000"}
    OUT.write("""DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL\tDrug-Resistance
COLOR\t#ff0000

LEGEND_TITLE\tDrug resistance
LEGEND_SHAPES\t1\t1\t1\t1
LEGEND_COLORS\t#80FF00\t#\t    #8000FF\t#FF0000
LEGEND_LABELS\tSensitive\tDrug-resisant\tMDR\tXDR

DATA
""")
    for s in samples:
        OUT.write("%s\t%s\n" % (s,dr_cols.get(results[s]["drtype"],"#000000")))
    OUT.close()

    drug_list = ['rifampicin', 'isoniazid', 'ethambutol', 'pyrazinamide', 'streptomycin', 'fluoroquinolones', 'aminoglycosides', 'kanamycin', 'amikacin', 'capreomycin', 'ethionamide', 'para-aminosalicylic_acid', 'clofazimine', 'linezolid', 'bedaquiline']
    OUT = open(prefix+".dr.indiv.itol.txt","w")
    dr_cols = {"Sensitive":"#80FF00","Drug-resistant":"#00FFFF","MDR":"#8000FF","XDR":"#FF0000"}
    legend_shapes = "\t".join(["2" for x in drug_list])
    legend_colours = "\t".join(["black" for x in drug_list])
    legend_labels = "\t".join(drug_list)
    OUT.write("""DATASET_BINARY
SEPARATOR TAB
DATASET_LABEL\tDrugs
COLOR\t#ff0000

SHOW_LABELS\t1
FIELD_SHAPES\t%s
FIELD_COLORS\t%s
FIELD_LABELS\t%s

DATA
""" % (legend_shapes,legend_colours,legend_labels))
    for s in samples:
        OUT.write("%s\t%s\n" % (s,"\t".join(["1" if d in dr_drugs[s] else "0" for d in drug_list])))
    OUT.close()
