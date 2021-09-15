import sys
import pathogenprofiler as pp
import json
from collections import defaultdict
from .utils import get_genome_positions_from_json_db, get_lt2drugs,rv2genes


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
        for d in x["drugs"]:
            drug = d["drug"]
            if float(x["freq"])<reporting_af:continue
            if drug not in results: results[drug] = []
            results[d["drug"]].append("%s %s (%.2f)" % (x["gene"],x["change"],float(x["freq"])))
            if drug not in annotation: annotation[drug] = {key:[] for key in columns}
            for key in columns:
                annotation[drug][key].append(x["drugs"][drug][key])
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

def select_most_relevant_csq(csqs):
    rank = ["frameshift_variant","start_lost","disruptive_inframe_deletion","disruptive_inframe_insertion","stop_gained","conservative_inframe_deletion","conservative_inframe_insertion","missense_variant","non_coding_transcript_exon_variant","upstream_gene_variant","stop_retained_variant","synonymous_variant"]
    ranked_csq = []
    for csq in csqs:
        ranked_csq.append([i for i,d in enumerate(rank) if d in csq["type"]][0])
    
    return csqs[ranked_csq.index(min(ranked_csq))]

def select_csq(dict_list):
    for d in dict_list:
        print(d)
        annotated_csq = []
        for csq in d["consequences"]:
            if "annotation" in csq:
                annotated_csq.append(csq)
        if len(annotated_csq)==0:
            csq = select_most_relevant_csq(d["consequences"])
            alternate_consequences = [json.dumps(x) for x in d["consequences"]]
            alternate_consequences.remove(json.dumps(csq))
            alternate_consequences = [json.loads(x) for x in alternate_consequences]
        elif len(annotated_csq)==1:
            csq = annotated_csq[0]
            alternate_consequences = []
        else:
            quit("ERROR! too many csqs")
        del d["consequences"]
        d.update(csq)
        d["alternate_consequences"] = alternate_consequences
    return dict_list

def dict_list_add_genes(dict_list,conf):
    rv2gene = {}
    for l in open(conf["bed"]):
        row = l.rstrip().split()
        rv2gene[row[3]] = row[4]
    for d in dict_list:
        print(d)
        d["locus_tag"] = d["gene_id"]
        d["gene"] = rv2gene[d["gene_id"]]
        del d["gene_id"]
        if "gene_name" in d:
            del d["gene_name"]
    return dict_list

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


def reformat_annotations(results,conf,reporting_af=0.1):
    #Chromosome      4998    Rv0005  -242
    lt2drugs = get_lt2drugs(conf["bed"])
    chr2gene_pos = {}
    for l in open(conf["ann"]):
        row = l.rstrip().split()
        chr2gene_pos[int(row[1])] = int(row[3])
    # for var in results["variants"]:
        # var["_internal_change"] = var["change"]
        # var["change"] = pp.reformat_mutations(var["change"],var["type"],var["locus_tag"],chr2gene_pos)
    resistant_drugs = set()
    results["dr_variants"] = []
    results["other_variants"] = []
    for var in results["variants"]:
        if "annotation" in var:
            tmp = var.copy()
            drvar = any([x["type"]=="drug" for x in var["annotation"]])
            phylovar = any([x["type"]=="phylogenetic" for x in var["annotation"]])
            if drvar:
                tmp["drugs"] = var["annotation"]
                del tmp["annotation"]
                if tmp["freq"]>=reporting_af:
                    for d in tmp["drugs"]:
                        resistant_drugs.add(d["drug"])
                results["dr_variants"].append(tmp)
            elif phylovar:
                var["lineage_variant"] = var["annotation"][0]["lineage"]
                del var["annotation"]
                results["other_variants"].append(var)
        else:
            var["gene_associated_drugs"] = lt2drugs[var["locus_tag"]]
            results["other_variants"].append(var)
    del results["variants"]

    FLQ_set = set(["moxifloxacin","levofloxacin","ciprofloxacin","ofloxacin"])
    SLI_set = set(["amikacin","capreomycin","kanamycin"])

    rif = "rifampicin" in resistant_drugs
    inh = "isoniazid" in resistant_drugs
    flq = len(FLQ_set.intersection(resistant_drugs)) > 0
    sli = len(SLI_set.intersection(resistant_drugs)) > 0

    if len(resistant_drugs)==0:
        drtype = "Sensitive"
    elif (rif and not inh) or (inh and not rif):
        drtype = "Pre-MDR"
    elif (rif and inh) and (not flq and not sli):
        drtype = "MDR"
    elif (rif and inh) and ( (flq and not sli) or (sli and not flq) ):
        drtype = "Pre-XDR"
    elif (rif and inh) and (flq and sli):
        drtype = "XDR"
    else:
        drtype = "Other"


    results["drtype"] = drtype
    return results

unlist = lambda t: [item for sublist in t for item in sublist]

def reformat_missing_genome_pos(results,conf):
    rv2gene = rv2genes(conf["bed"])
    dr_associated_genome_pos = get_genome_positions_from_json_db(conf["json_db"],conf["ann"],conf["gff"])
    genome_pos2gene_pos = {}
    genome_pos2gene = {}
    for l in open(conf["ann"]):
        row = l.strip().split()
        genome_pos2gene_pos[int(row[1])] = int(row[3])
        genome_pos2gene[int(row[1])] = row[2]

    tmp_results = defaultdict(lambda: defaultdict(list))
    for genome_pos in results:
        # This is a bit dangerous as it assumes all coding genes are labled with "Rv"
        gene = genome_pos2gene[genome_pos]
        gene_pos = genome_pos2gene_pos[genome_pos]
        coding = True if gene[:2]=="Rv" else False
        if coding:
            if gene_pos>=0:
                codon = ((gene_pos-1)//3) + 1
                tmp_results[gene][codon].append(genome_pos)
            else:
                tmp_results[gene][gene_pos].append(genome_pos)
        else:
            tmp_results[gene][gene_pos].append(genome_pos)
    new_results = []
    for gene in tmp_results:
        for pos in tmp_results[gene]:
            genome_positions = tmp_results[gene][pos]
            dr_position = list(set(unlist([unlist([y[2] for y in dr_associated_genome_pos[x]]) for x in genome_positions if x in dr_associated_genome_pos])))

            new_results.append({"locus_tag":gene, "gene": rv2gene[gene], "genome_positions": genome_positions , "position":pos, "position_type":"codon" if (gene[:2]=="Rv" and pos>=0) else "gene", "drug_resistance_position": dr_position})
    return new_results



def reformat(results,conf,reporting_af,mutation_metadata=False):
    results["variants"] = [x for x in results["variants"] if len(x["consequences"])>0]
    results["variants"] = select_csq(results["variants"])
    results["variants"] = dict_list_add_genes(results["variants"],conf)
    if "gene_coverage" in results["qc"]:
        results["qc"]["gene_coverage"] = dict_list_add_genes(results["qc"]["gene_coverage"],conf)
        results["qc"]["missing_positions"] = reformat_missing_genome_pos(results["qc"]["missing_positions"],conf)
    results = barcode2lineage(results)
    results = reformat_annotations(results,conf,reporting_af)
    results["db_version"] = json.load(open(conf["version"]))
    if mutation_metadata:
        pass
        # results = add_mutation_metadata(results)
    return results
