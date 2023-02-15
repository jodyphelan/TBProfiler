import sys
import argparse
from docxtpl import DocxTemplate
import json
from collections import defaultdict
import datetime


def sanitize(d):
    d = d.replace("-","_")
    return d


def write_docx(json_results,conf,outfile,template_file,reporting_af = 0.0):
    data = json_results
    drug_variants = defaultdict(list)
    confidence = defaultdict(list)
    for var in data['dr_variants']:
        for d in var['drugs']:
            if d['type']=='drug':
                if "who confidence" not in d: continue
                drug_variants[d['drug']].append(f"{var['gene']}_{var['change']}")
                confidence[d['drug']].append(d['who confidence'])

    
    time = datetime.datetime.strptime(data['timestamp'], "%d-%m-%Y %H:%M:%S")

    variables = {
        'date':time.strftime("%d %b %Y"),
        'sublineage': data['sublin'],
        'resistant_drugs': ", ".join([d.capitalize() for d in conf['drugs'] if d in drug_variants]),
        'version': data['tbprofiler_version'],
        'sensitive': True if data['drtype'] == "Sensitive" else False,
        'mdr': True if data['drtype'] in ("MDR-TB","Pre-XDR-TB") else False,
        'xdr': True if data['drtype'] in "XDR-TB" else False,
        'resistant': True if data['drtype'] in ("Other","RR-TB","HR-TB") else False
    }

    for d in conf['drugs']:
        variables[sanitize(d)+"_variants"] = ", ".join(drug_variants[d]) if len(drug_variants[d])>0 else "Not found"
        variables[sanitize(d)+"_confidence"] = ", ".join(confidence[d]) if len(drug_variants[d])>0 else "-"
        variables[sanitize(d)+"_interpretation"] = "Resistance" if len(drug_variants[d])>0 else "-"


    doc = DocxTemplate(template_file)
    doc.render(variables)
    doc.save(outfile)