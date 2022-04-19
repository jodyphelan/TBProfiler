import sys
import pathogenprofiler as pp
import json
from collections import defaultdict
from .utils import get_genome_positions_from_json_db, get_lt2drugs,rv2genes
from .xdb import *




def get_main_lineage(lineage_dict_list,max_node_skip=1):
    def collapse_paths(paths):
        filtered_paths = []
        for p in sorted(paths,reverse=True):
            if p=="lineageBOV_AFRI": continue
            path_stored = any([p in x for x in filtered_paths])
            if not path_stored:
                filtered_paths.append(p)
        return filtered_paths

    def derive_path(x):
        return [".".join(x.split(".")[:i])for i in range(1,len(x.split(".")))] + [x]

    lin_freqs = {}
    pool = []
    for l in lineage_dict_list:
        pool.append(l["lin"].replace("M.","M_"))
        lin_freqs[l["lin"].replace("M.","M_")] = float(l["frac"])
    routes = [";".join(derive_path(x)) for x in pool]
    paths = collapse_paths(routes)
    path_mean_freq = {}
    for path in paths:
        nodes = tuple(path.split(";"))
        nodes_skipped = sum([n not in pool for n in nodes])
        if nodes_skipped>max_node_skip: continue
        freqs = [lin_freqs[n] for n in nodes if n in lin_freqs]
        path_mean_freq[nodes] = sum(freqs)/len(freqs)
    main_lin = ";".join(sorted(list(set([x[0] for x in path_mean_freq])))).replace("_",".")
    sublin = ";".join([x[-1] for x in path_mean_freq]).replace("_",".")
    return (main_lin,sublin)

def barcode2lineage(results,max_node_skip=1):
    results["lineage"] = []
    for d in results["barcode"]:
        results["lineage"].append({"lin":d["annotation"],"family":d["info"][0],"spoligotype":d["info"][1],"rd":d["info"][2],"frac":d["freq"]})
    del results["barcode"]
    results["lineage"] = sorted(results["lineage"],key= lambda x:len(x["lin"]))
    main_lin,sublin = get_main_lineage(results["lineage"])
    results["main_lin"] = main_lin
    results["sublin"] = sublin
    return results




def add_drtypes(results,reporting_af=0.1):
    resistant_drugs = set()
    for var in results["dr_variants"]:
        if var["freq"]>=reporting_af:
            for d in var["drugs"]:
                resistant_drugs.add(d["drug"])

    FLQ_set = set(["levofloxacin","moxifloxacin","ciprofloxacin","ofloxacin"])
    groupA_set = set(["bedaquiline","linezolid"])

    rif = "rifampicin" in resistant_drugs
    inh = "isoniazid" in resistant_drugs
    flq = len(FLQ_set.intersection(resistant_drugs)) > 0
    gpa = len(groupA_set.intersection(resistant_drugs)) > 0

    if len(resistant_drugs)==0:
        drtype = "Sensitive"
    elif (rif and not inh) and not flq:
        drtype = "RR-TB"
    elif (inh and not rif):
        drtype = "HR-TB"
    elif (rif and inh) and not flq:
        drtype = "MDR-TB"
    elif rif and (flq and not gpa):
        drtype = "Pre-XDR-TB"
    elif rif and (flq and gpa):
        drtype = "XDR-TB"
    else:
        drtype = "Other"


    results["drtype"] = drtype
    return results

unlist = lambda t: [item for sublist in t for item in sublist]



def suspect_profiling(results):
    new_vars = []
    for var in results["other_variants"]:
        if var["type"]!="missense_variant": continue
        pred = None
        if var["gene"]=="atpE":
            pred = get_biosig_bdq_prediction(var["change"])
        if var["gene"]=="pncA":
            pred = get_biosig_pza_prediction(var["change"])
        if pred:
            if "annotation" in var:
                var["annotation"].append(pred)
            else:
                var["annotation"] = [pred]
            var["drugs"] = [{
                "type":"drug",
                "drug":"pyrazinamide" if var["gene"]=="pncA" else "bedquiline",
                "confers": "resistance",
                "evidence": "suspect-PZA" if var["gene"]=="pncA" else "suspect-BDQ"
            }]
            new_vars.append(var)
    for v in new_vars:
        results["dr_variants"].append(v)
        results["other_variants"].remove(v)
    return results

    
def reformat(results,conf,reporting_af,mutation_metadata=False,use_suspect=False):
    results["variants"] = [x for x in results["variants"] if len(x["consequences"])>0]
    results["variants"] = pp.select_csq(results["variants"])
    results["variants"] = pp.dict_list_add_genes(results["variants"],conf)
    if "gene_coverage" in results["qc"]:
        results["qc"]["gene_coverage"] = pp.dict_list_add_genes(results["qc"]["gene_coverage"],conf)
    if "missing_positions" in results["qc"]:
        results["qc"]["missing_positions"] = pp.reformat_missing_genome_pos(results["qc"]["missing_positions"],conf)
    if "barcode" in results:
        results = barcode2lineage(results)
    results = pp.reformat_annotations(results,conf)
    results = add_drtypes(results,reporting_af)
    results["db_version"] = json.load(open(conf["version"]))
    if mutation_metadata:
        pass
        # results = add_mutation_metadata(results)
    if use_suspect:
        results = suspect_profiling(results)
    return results
