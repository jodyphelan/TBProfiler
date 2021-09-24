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
    return g.json()

def get_biosig_bdq_prediction(raw_change):
    change = aa_long2short(raw_change)
    g = requests.get(f"http://biosig.unimelb.edu.au/suspect_bdq/api/prediction_single?mutation={change}")
    return g.json()