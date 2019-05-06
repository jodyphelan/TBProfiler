import pathogenprofiler as pp
import time
from .reformat import *

def dict_list2html(l,columns = None, mappings = None):
    headings = list(l[0].keys()) if not columns else columns
    rows = []
    header = "<tr>%s</tr>" % "".join(["<td>%s</td>" % mappings[x] if (mappings!=None and x in mappings) else "<td>%s</td>" % x.title() for x in headings])
    for row in l:
        r = "<tr>%s</tr>" % "".join(["<td>%s</td>" % "%.3f" % row[x] if isinstance(row[x],float) else "<td>%s</td>" % str(row[x]).replace("_", " ") for x in headings])
        rows.append(r)
    str_rows = "\n".join(rows)
    return "<table cellpadding=0 cellspacing=0 width=100%%>\n%s\n%s\n</table>" % (header,str_rows)

def load_html(html_strings):
	return r"""
<p>
<p>
<a href="input/results/%(id)s.results.csv">Download results <img src="/img/download_icon.jpeg" height="20px" width="20px"></img></a>
<div><b>Summary</b></div>
<div><b>ID:</b> %(id)s</div>
<div><b>Date:</b> %(date)s</div>
<div><b>Strain:</b> %(strain)s</div>
<div><b>Drug-resistance:</b> %(drtype)s</div>
</p>
<p>
%(dr_report)s
</p>
<p>
%(lineage_report)s
</p>
<p>
<div>Mutations in candidate genes:</div>
%(other_var_report)s
</p>

<p>
<div>Analysis pipeline specifications<div>
<p><b>TBProfiler Version:</b> %(version)s</p>
%(pipeline)s
</p>
""" % html_strings

def write_html(json_results,conf,outfile,columns = None,drug_order = None):
	json_results = get_summary(json_results,conf,columns = columns, drug_order=drug_order)
	html_strings = {}
	html_strings["id"] = json_results["id"]
	html_strings["date"] = time.ctime()
	html_strings["strain"] = json_results["sublin"]
	html_strings["drtype"] = json_results["drtype"]
	html_strings["dr_report"] = dict_list2html(json_results["drug_table"],["Drug","Genotypic Resistance","Mutations"]+columns,{"Drug":"Drug<sup>1</sup>","Genotypic Resistance":"Resistance","Mutations":"Supporting Mutations (frequency)"})
	html_strings["lineage_report"] = dict_list2html(json_results["lineage"],["lin","family","spoligotype","rd"],{"lin":"Lineage<sup>2</sup>","frac":"Estimated fraction","family":"Family","spoligotype":"Main Spoligotype","rd":"RDS"})
	html_strings["other_var_report"] = dict_list2html(json_results["other_variants"],["gene","genome_pos","change","freq"],{"gene":"Gene","genome_pos":"Chromosome Position","change":"Mutation","freq":"Estimated fraction"})
	html_strings["pipeline"] = dict_list2html(json_results["pipline_table"],["Analysis","Program"])
	html_strings["version"] = json_results["tbprofiler_version"]
	o = open(outfile,"w")
	pp.log("Writing results to %s" % outfile)
	o.write(load_html(html_strings))
	o.close()
