import pathogenprofiler as pp
import time
def dict_list2html(l,columns = None, mappings = {}):
    headings = list(l[0].keys()) if not columns else columns
    rows = []
    header = "<tr>%s</tr>" % "".join(["<td>%s</td>" % mappings[x] if x in mappings else "<td>%s</td>" % x.title() for x in headings])
    for row in l:
        r = "<tr>%s</tr>" % "".join(["<td>%s</td>" % "%.3f" % row[x] if type(row[x])==float else "<td>%s</td>" % str(row[x]).replace("_", " ") for x in headings])
        rows.append(r)
    str_rows = "\n".join(rows)
    return "<table cellpadding=0 cellspacing=0 width=100%%>\n%s\n%s\n</table>" % (header,str_rows)

def load_html(html_strings):
	return r"""
<p>
<p>
<a href="input/results/%(id)s.results.csv">Download results <img src="/img/download_icon.jpeg"></img></a>
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

def write_html(json_results,conf,outfile,columns = [],drug_order = None):
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
		if d not in results: results[d] = list()
		results[d].append("%s %s (%.2f)" % (x["gene"],x["change"].replace(">","&gt;"),x["freq"]))
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
	html_strings = {}
	html_strings["id"] = json_results["id"]
	html_strings["date"] = time.ctime()
	html_strings["strain"] = json_results["sublin"]
	html_strings["drtype"] = json_results["drtype"]
	html_strings["dr_report"] = dict_list2html(drug_table,["Drug","Genotypic Resistance","Mutations"]+columns,{"Drug":"Drug<sup>1</sup>","Genotypic Resistance":"Resistance","Mutations":"Supporting Mutations (frequency)"})
	html_strings["lineage_report"] = dict_list2html(json_results["lineage"],["lin","family","spoligotype","rd"],{"lin":"Lineage<sup>2</sup>","frac":"Estimated fraction","family":"Family","spoligotype":"Main Spoligotype","rd":"RDS"})
	html_strings["other_var_report"] = dict_list2html(json_results["other_variants"],["gene","genome_pos","change","freq"],{"gene":"Gene","genome_pos":"Chromosome Position","change":"Mutation","freq":"Estimated fraction"})
	html_strings["pipeline"] = dict_list2html(pipeline_tbl,["Analysis","Program"])
	html_strings["version"] = json_results["tbprofiler_version"]
	o = open(outfile,"w")
	pp.log("Writing results to %s" % outfile)
	o.write(load_html(html_strings))
	o.close()
