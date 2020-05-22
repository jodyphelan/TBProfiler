import re
import pathogenprofiler as pp
import json

def get_summary(json_results,conf,columns = None,drug_order = None,reporting_af=0.0):
    if not columns:
        columns=[]
    drugs = set()
    for l in open(conf["bed"]):
        arr = l.rstrip().split()
        for d in arr[5].split(","):
            drugs.add(d)
    if drug_order:
        drugs = drug_order
    drug_table = []
    results = {}
    annotation = {}
    for key in columns:
        if key not in json_results["dr_variants"][0]: pp.log("%s not found in variant annotation, is this a valid column in the database CSV file? Exiting!" % key,True)
    for x in json_results["dr_variants"]:
        d = x["drug"]
        if float(x["freq"])<reporting_af:continue
        if d not in results: results[d] = list()
        results[d].append("%s %s (%.2f)" % (x["gene"],x["change"],x["freq"]))
        if d not in annotation: annotation[d] = {key:[] for key in columns}
        for key in columns:
            annotation[d][key].append(x[key])
    for d in drugs:
        if d in results:
            results[d] = ", ".join(results[d]) if len(results[d])>0 else ""
            r = "R" if len(results[d])>0 else ""
            for key in columns:
                annotation[d][key] = ", ".join(annotation[d][key]) if len(annotation[d][key])>0 else ""
        else:
            results[d] = ""
            r = ""
        dictline = {"Drug":d.capitalize(),"Genotypic Resistance":r,"Mutations":results[d]}
        for key in columns:
            dictline[key] = annotation[d][key] if d in annotation else ""
        drug_table.append(dictline)
    pipeline_tbl = [{"Analysis":"Mapping","Program":json_results["pipeline"]["mapper"]},{"Analysis":"Variant Calling","Program":json_results["pipeline"]["variant_caller"]}]
    new_json = json_results.copy()
    new_json["drug_table"] = drug_table
    new_json["pipline_table"] = pipeline_tbl
    return new_json

def dict_list_add_genes(dict_list,conf):
    rv2gene = {}
    for l in open(conf["bed"]):
        row = l.rstrip().split()
        rv2gene[row[3]] = row[4]
    for d in dict_list:
        d["locus_tag"] = d["gene_id"]
        d["gene"] = rv2gene[d["gene_id"]]
        del d["gene_id"]
    return dict_list

def add_genes(results,conf):

    rv2gene = {}
    for l in open(conf["bed"]):
        row = l.rstrip().split()
        rv2gene[row[3]] = row[4]
    for d in results["variants"]:
        d["locus_tag"] = d["gene_id"]
        d["gene"] = rv2gene[d["gene_id"]]
        del d["gene_id"]
    return results

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
        pool.append(l["lin"])
        lin_freqs[l["lin"]] = float(l["frac"])
    routes = [";".join(derive_path(x)) for x in pool]
    paths = collapse_paths(routes)
    path_mean_freq = {}
    for path in paths:
        nodes = tuple(path.split(";"))
        nodes_skipped = sum([n not in pool for n in nodes])
        if nodes_skipped>max_node_skip: continue
        freqs = [lin_freqs[n] for n in nodes if n in lin_freqs]
        path_mean_freq[nodes] = sum(freqs)/len(freqs)

    main_lin = ";".join(sorted(list(set([x[0] for x in path_mean_freq]))))
    sublin = ";".join([x[-1] for x in path_mean_freq])
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


def reformat_annotations(results,conf,reporting_af=0.1):
    #Chromosome      4998    Rv0005  -242
    chr2gene_pos = {}
    for l in open(conf["ann"]):
        row = l.rstrip().split()
        chr2gene_pos[int(row[1])] = int(row[3])
    for var in results["variants"]:
        var["_internal_change"] = var["change"]
        var["change"] = pp.reformat_mutations(var["change"],var["type"],var["locus_tag"],chr2gene_pos)
    results["dr_variants"] = []
    for d in [x for x in results["variants"] if "annotation" in x]:
        for drug in d["annotation"]["drugs"]:
            tmp = d.copy()
            tmp["drug"] = drug
            for key in d["annotation"]["drugs"][drug]:
                tmp[key] = d["annotation"]["drugs"][drug][key]
            del tmp["annotation"]
            results["dr_variants"].append(tmp)
    results["other_variants"] = [x for x in results["variants"] if "annotation" not in x]
    del results["variants"]
    dr_drugs = [x["drug"] for x in results["dr_variants"] if x["freq"]>=reporting_af]
    MDR = "R" if ("isoniazid" in dr_drugs and "rifampicin" in dr_drugs) else ""
    XDR = "R" if MDR=="R" and ( "amikacin" in dr_drugs or "kanamycin" in dr_drugs or "capreomycin" in dr_drugs ) and ( "fluoroquinolones" in dr_drugs) else ""
    drtype = "Sensitive"
    if XDR=="R":
        drtype="XDR"
    elif MDR=="R":
        drtype="MDR"
    elif len(dr_drugs)>0:
        drtype="Drug-resistant"
    results["XDR"] = XDR
    results["MDR"] = MDR
    results["drtype"] = drtype
    return results

def reformat(results,conf,reporting_af):
    results["variants"] = dict_list_add_genes(results["variants"],conf)
    if "gene_coverage" in results["qc"]:
        results["qc"]["gene_coverage"] = dict_list_add_genes(results["qc"]["gene_coverage"],conf)
    results = barcode2lineage(results)
    results = reformat_annotations(results,conf,reporting_af)
    results["db_version"] = json.load(open(conf["version"]))
    return results
