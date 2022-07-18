import pathogenprofiler as pp
import requests
import re


def aa_long2short(mut):
    aconv = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C",
    "Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I",
    "Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P",
    "Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V",
    "Stop":"*", "-":"-","*":"*"
    }
    r = re.search("p.([A-Z][a-z][a-z])([0-9]+)([A-Za-z\*]+)",mut)
    return f"{aconv[r.group(1)]}{r.group(2)}{aconv[r.group(3)]}"

def get_biosig_pza_prediction(raw_change):
    change = aa_long2short(raw_change)
    g = requests.get(f"http://biosig.unimelb.edu.au/suspect_pza/api/prediction_single?mutation={change}")
    data = g.json()
    data["prediction"] = data["suspect_pza_prediction"]
    return data

def get_biosig_bdq_prediction(raw_change):
    change = aa_long2short(raw_change)
    g = requests.get(f"http://biosig.unimelb.edu.au/suspect_bdq/api/prediction_single?mutation={change}")
    data = g.json()
    data["prediction"] = data["suspect_bdq_prediction"]
    return data

def suspect_profiling(results):
    new_vars = []
    for var in results["other_variants"]:
        if var["type"]!="missense_variant": continue
        pred = None
        if var["gene"]=="atpE":
            pp.infolog(f"Profiling {var['gene']} {var['change']} with suspect-bdq")
            pred = get_biosig_bdq_prediction(var["change"])
        if var["gene"]=="pncA":
            pp.infolog(f"Profiling {var['gene']} {var['change']} with suspect-pza")
            pred = get_biosig_pza_prediction(var["change"])
        if pred:
            if "annotation" in var:
                var["annotation"].append(pred)
            else:
                var["annotation"] = [pred]
            if pred["prediction"]=="Resistant":
                var["drugs"] = [{
                    "type":"drug",
                    "drug":"pyrazinamide" if var["gene"]=="pncA" else "bedaquiline",
                    "confers": "resistance",
                    "evidence": "suspect-PZA" if var["gene"]=="pncA" else "suspect-BDQ"
                }]
                new_vars.append(var)
    for v in new_vars:
        results["dr_variants"].append(v)
        results["other_variants"].remove(v)
    return results