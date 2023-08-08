import json
import os
from collections import defaultdict
from tqdm import tqdm
from .utils import get_lt2drugs
import logging
import csv 

def get_common_fields(rows):
    cols = set(rows[0].keys())
    orders = {}
    for r in rows:
        cols = cols.intersection(set(r.keys()))
        for i,c in enumerate(r):
            orders[c] = i
    return sorted(cols,key=lambda x:orders[x])

def get_field_values(rows,cols):
    new_rows = []
    for r in rows:
        new_rows.append({c:r.get(c,None) for c in cols})
    return new_rows


def collate_results(prefix,conf,result_dirs=["./results"],sample_file=None,full_results=True,full_variant_results=True,mark_missing=False,sep="\t"):
    for d in result_dirs:
        if not os.path.isdir(d):
            errlog("\nERROR: Can't find directory %s\n" % d )
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


    samples = {}
    for d in result_dirs:
        for f in os.listdir(f"{d}/"):
            if f.endswith(".results.json"):
                s = f[:-13]
                samples[s] = f'{d}/{s}.results.json'

    if sample_file:
        samples = {s:samples[s] for s in [l.strip() for l in open(sample_file)]}

    results = defaultdict(dict)
    result_rows = []
    dr_variants = defaultdict(lambda:defaultdict(dict))
    dr_variants_set = set()
    dr_drugs = {}
    sample_dr_mutations_set = defaultdict(set)
    sample_other_mutations_set = defaultdict(set)
    for s in samples:
        for d in drug_list:
            results[s][d] = set()
    lt2drugs = get_lt2drugs(conf["bed"])

    tmp_edges = set()
    for s in tqdm(samples):
        res = {"sample":s}
        dr_drugs[s] = set()
        temp = json.load(open(samples[s]))

        missing_drugs = set()
        if "missing_positions" in temp["qc"]:
            missing_genes = [x["locus_tag"] for x in temp["qc"]["missing_positions"] if "drugs" in x]
            for locus_tag in missing_genes:
                missing_drugs = missing_drugs.union(set(lt2drugs[locus_tag]))


        for x in temp["dr_variants"]:
            dr_variants[x["gene"]][x["change"]][s] = x["freq"]
            sample_dr_mutations_set[s].add((x["gene"],x["change"]))
            dr_variants_set.add((x["gene"],x["change"]))
            for d in x["drugs"]:
                dr_drugs[s].add(d["drug"])
                results[s][d["drug"]].add("%s_%s" % (x["gene"],x["change"]) if full_results else "R")
        for x in temp["other_variants"]:
            sample_other_mutations_set[s].add((x["gene"],x["change"]))
        
        res["main_lineage"] = results[s]["main_lin"] = temp["main_lin"]
        res["sub_lineage"] = results[s]["sublin"] = temp["sublin"]
        if "spoligotype" in temp:
            res["spoligotype"] = temp["spoligotype"]["octal"]
        res["DR_type"] = results[s]["drtype"] = temp["drtype"]
        res["region_median_depth"] = results[s]["region_median_depth"] = temp["qc"].get("region_median_depth","NA")
        res["pct_reads_mapped"] = temp["qc"].get("pct_reads_mapped","NA")
        res["num_reads_mapped"] = temp["qc"].get("num_reads_mapped","NA")
        res["num_dr_variants"] = len(sample_dr_mutations_set[s])
        res["num_other_variants"] = len(sample_other_mutations_set[s])
        for d in drug_list:
            results[s][d] = ", ".join(sorted(results[s][d])) if len(results[s][d])>0 else "-"
            if mark_missing and d in missing_drugs:
                res[d] = "*%s" % results[s][d]
            else:
                res[d] = results[s][d]
        result_rows.append(res)


        if "close_samples" in temp:
            for d in temp['close_samples']:
                sorted_samples = sorted([s,d['sample']])
                tmp_edges.add((sorted_samples[0],sorted_samples[1],d['distance']))

    with open(prefix+".txt","w") as OUT:
        fields = get_common_fields(result_rows)
        writer = csv.DictWriter(OUT,fieldnames=fields,delimiter=sep)
        writer.writeheader()
        writer.writerows(get_field_values(result_rows,fields))

    if full_variant_results:

        all_vars = conf["json_db"]
        lt2gene = {}
        for l in open(conf["bed"]):
            #Chromosome      5240    7267    Rv0005  gyrB    FLUOROQUINOLONES
            row = l.rstrip().split()
            lt2gene[row[3]] = row[4] if row[4]!="." else row[3]
        for gene in all_vars:
            for mutation in all_vars[gene]:
                dr_variants_set.add((lt2gene[gene],mutation))
        VAR = open(prefix+".variants.txt","w")
        VAR.write("sample\t%s\n" % ("\t".join(["%s_%s" % (g,c) for g,c in sorted(dr_variants_set,key=lambda x: x[0])])))
        for s in samples:
            VAR.write("%s\t%s\n" % (s,"\t".join(["%.3f" % dr_variants[gene][change][s] if gene in dr_variants and change in dr_variants[gene] and s in dr_variants[gene][change] else "0" for gene,change in sorted(dr_variants_set,key=lambda x: x[0])])))
        VAR.close()
    
    # OUT = 
    # OUT.write("sample\tmain_lineage\tsub_lineage\tDR_type\tpct_reads_mapped\tnum_reads_mapped\tmedian_coverage\tnum_dr_variants\tnum_other_variants\t%s" % "\t".join(drug_list)+"\n")
    # for s in samples:
        
        # OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(s,results[s]["main_lin"],results[s]["sublin"],results[s]["drtype"],results[s]["pct_reads_mapped"],results[s]["num_reads_mapped"],results[s]["median_coverage"],results[s]["num_dr_variants"],results[s]["num_other_variants"],"\t".join([results[s][x] for x in drug_list])))
    # OUT.close()
    json.dump(dict(results),open(prefix+".json","w"))
    
    
    
    lineage_cols = {"lineage1":"#104577","lineage2":"#ab2323","lineage3":"#18a68c","lineage4":"#f68e51","lineage5":"#7cb5d2","lineage6":"#fde05e","lineage7":"#bc94b7","lineage8":"#ccc9e7","lineage9":"#bd9391","Animal strains":"#f8e0c8","Other":"#000000"}
    lineage_aggregation = {"M.caprae":"Animal strains","M.bovis":"Animal strains","M.orygis":"Animal strains"}
    lineages_present = set([lineage_aggregation.get(results[s]["main_lin"],results[s]["main_lin"]) if ";" not in results[s]["main_lin"] else "Other" for s in samples])
    lineage_cols = {key:val for key,val in lineage_cols.items() if key in lineages_present}
    OUT = open(prefix+".lineage.itol.txt","w")
    OUT.write("""DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL\tLineage
COLOR\t#ff0000

LEGEND_TITLE\tLineage
LEGEND_SHAPES\t%s
LEGEND_LABELS\t%s
LEGEND_COLORS\t%s

DATA
""" % ("\t".join(["1" for _ in lineage_cols]),"\t".join(list([x.title() for x in lineage_cols.keys()])), "\t".join(list(lineage_cols.values()))))
    for s in samples:
        OUT.write("%s\t%s\n" % (s,lineage_cols.get(lineage_aggregation.get(results[s]["main_lin"],results[s]["main_lin"]),"#000000")))
    OUT.close()

    OUT = open(prefix+".dr.itol.txt","w")
    dr_cols = {"Sensitive":"#28a745","RR-TB":"#007bff","HR-TB":"#E0ACD5","MDR-TB":"#ffc107","Pre-XDR-TB":"#dc3545","XDR-TB":"#343a40","Other":"#f8f9fa"}
    drtypes_present = set([results[s]["drtype"] for s in samples])
    dr_cols = {key:val for key,val in dr_cols.items() if key in drtypes_present}
    OUT.write("""DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL\tDrug-Resistance
COLOR\t#ff0000

LEGEND_TITLE\tDrug resistance

LEGEND_SHAPES\t%s
LEGEND_LABELS\t%s
LEGEND_COLORS\t%s


DATA
""" % ("\t".join(["1" for _ in dr_cols]),"\t".join(list(dr_cols.keys())), "\t".join(list(dr_cols.values()))))
    for s in samples:
        OUT.write("%s\t%s\n" % (s,dr_cols.get(results[s]["drtype"],"#000000")))
    OUT.close()

    drug_list = ['rifampicin', 'isoniazid', 'ethambutol', 'pyrazinamide', 'streptomycin', 'fluoroquinolones', 'aminoglycosides', 'kanamycin', 'amikacin', 'capreomycin', 'ethionamide', 'para-aminosalicylic_acid', 'clofazimine', 'linezolid', 'bedaquiline', 'delamanid']
    OUT = open(prefix+".dr.indiv.itol.txt","w")
    # dr_cols = {"Sensitive":"#80FF00","Drug-resistant":"#00FFFF","MDR":"#8000FF","XDR":"#FF0000"}
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

    if len(tmp_edges)>0:
        tmp_nodes = set()
        edges = []

        for i,e in enumerate(tmp_edges):
            if e[0] in results and e[1] in results:
                tmp_nodes.add(e[0])
                tmp_nodes.add(e[1])
                edges.append({"source":e[0],"target":e[1],"properties":{"distance":e[2]}})
        nodes = []
        for n in tmp_nodes:
            nodes.append({"id":n,"properties":{"drtype":results[n]["drtype"],"lineage":results[n]["main_lin"],"median_depth":results[n]["region_median_depth"]}})
        json.dump({"nodes":nodes,"edges":edges},open(prefix+".transmission_graph.json","w"))

        with open(prefix+".distance_matrix.txt","w") as MAT:
            MAT.write("samples\t%s\n" % "\t".join(samples))
            transformed_edges = {}
            for e in tmp_edges:
                transformed_edges[(e[0],e[1])] = e[2]
            for si in samples:
                row = [si]
                for sj in samples:
                    ss = tuple(sorted([si,sj]))
                    if si==sj:
                        row.append("0")
                    else:
                        row.append(str(transformed_edges[ss]) if ss in transformed_edges else "NA")
                MAT.write("%s\n" % "\t".join(row))