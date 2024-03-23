import sys
from docxtpl import DocxTemplate
from collections import defaultdict


def sanitize(d):
    d = d.replace("-","_")
    return d


def write_docx(result,conf,outfile,template_file = None):
    if template_file is None:
        template_file = sys.prefix+"/share/tbprofiler/default_template.docx"
    data = result.model_dump()
    drug_variants = defaultdict(list)
    confidence = defaultdict(list)
    for var in data['dr_variants']:
        for d in var['drugs']:
            if d['type']=='drug_resistance':
                drug_variants[d['drug']].append(f"{var['gene_name']}_{var['change']}")
                confidence[d['drug']].append(d['confidence'])



    variables = {
        'date':data['timestamp'].strftime("%d %b %Y"),
        'sublineage': data['sub_lineage'],
        'resistant_drugs': ", ".join([d.capitalize() for d in conf['drugs'] if d in drug_variants]),
        'version': data['pipeline']['software_version'],
        'sensitive': True if data['drtype'] == "Sensitive" else False,
        'mdr': True if data['drtype'] in ("MDR-TB","Pre-XDR-TB") else False,
        'xdr': True if data['drtype'] in "XDR-TB" else False,
        'resistant': True if data['drtype'] in ("Other","RR-TB","HR-TB") else False,
        'drtype': data['drtype'],
        'd': data,
    }




    for d in conf['drugs']:
        variables[sanitize(d)+"_variants"] = ", ".join(drug_variants[d]) if len(drug_variants[d])>0 else "Not found"
        variables[sanitize(d)+"_confidence"] = ", ".join(confidence[d]) if len(drug_variants[d])>0 else "-"
        variables[sanitize(d)+"_interpretation"] = "Resistance" if len(drug_variants[d])>0 else "-"


    doc = DocxTemplate(template_file)
    doc.render(variables)
    doc.save(outfile)