import json
from collections import defaultdict

class profiling_results:
	params = {}
	drugs = set()
	samples = []
	def __init__(self,conf_file,samples_file,prefix,stor_dir,full_results):
		tmp = json.load(open(conf_file))
		for x in tmp:
			self.params[x] = tmp[x]

		for l in open(self.params["dr_bed_file"]):
			arr = l.rstrip().split()
			for d in arr[5].split(";"):
				self.drugs.add(d)
		for s in open(samples_file):
			self.samples.append(s.rstrip())
		self.prefix = prefix
		self.stor_dir = stor_dir
		self.full_results = full_results

	def results2tab(self):
		results = defaultdict(dict)
		linresults = defaultdict(dict)
		for s in self.samples:
			for d in self.drugs:
				results[s][d] = set()
		for s in self.samples:
			temp = json.load(open("%s/results/%s.results.json" % (self.stor_dir,s)))
			for x in temp["small_variants_dr"]:
				for d in x["drug"].split(";"):
					results[s][d].add("%s_%s" % (x["gene"],x["change"]) if self.full_results else "R")
			for x in temp["del"]:
				for d in x["drug"].split(";"):
					results[s][d].add("large_deletion_%s" % x["gene"] if self.full_results else "R")
			for d in self.drugs:
				results[s][d] = ", ".join(results[s][d]) if len(results[s][d])>0 else "-"
			linresults[s]["main"] = sorted([x["lin"] for x in temp["lineage"]])[0] if len(temp["lineage"])>0 else "-"
			linresults[s]["sublin"] = sorted([x["lin"] for x in temp["lineage"]])[-1] if len(temp["lineage"])>0 else "-"
			dr_drugs = [x["drug"] for x in temp["small_variants_dr"]]

			MDR = "R" if ("ISONIAZID" in dr_drugs and "RIFAMPICIN" in dr_drugs) else "-"
			XDR = "R" if MDR=="R" and ( "AMIKACIN" in dr_drugs or "KANAMYCIN" in dr_drugs or "CAPREOMYCIN" in dr_drugs ) and ( "FLUOROQUINOLONES" in dr_drugs) else "-"
			drtype = "Sensitive"
			if XDR=="R":
				drtype="XDR"
			elif MDR=="R":
				drtype="MDR"
			elif len(dr_drugs)>0:
				drtype="Drug-resistant"
			results[s]["XDR"] = XDR
			results[s]["MDR"] = MDR
			results[s]["drtype"] = drtype

		o = open(self.prefix+".txt","w")
		o.write("sample\tmain_lineage\tsub_lineage\tDR_type\tMDR\tXDR\t%s" % "\t".join(self.drugs)+"\n")
		for s in self.samples:
			o.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(s,linresults[s]["main"],linresults[s]["sublin"],results[s]["drtype"],results[s]["MDR"],results[s]["XDR"],"\t".join([results[s][x] for x in self.drugs])))
		o.close()

		lineage_cols = {"lineage1":"#104577","lineage2":"#ab2323","lineage3":"#18a68c","lineage4":"#f68e51","lineage5":"#7cb5d2","lineage6":"#fde05e","lineage7":"#bc94b7","lineageBOV":"#f8e0c8","lineageOther":"#000000"}
		o = open(self.prefix+".lineage.itol.txt","w")
		o.write("""DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL	Lineage
COLOR	#ff0000

LEGEND_TITLE	Lineage
LEGEND_SHAPES	1	1	1	1	1	1	1	1	1
LEGEND_COLORS	%(lineage1)s	%(lineage2)s	%(lineage3)s	%(lineage4)s	%(lineage5)s	%(lineage6)s	%(lineage7)s	%(lineageBOV)s	%(lineageOther)s
LEGEND_LABELS	Lineage1	Lineage2	Lineage3	Lineage4	Lineage5	Lineage6	Lineage7	Bovis	Other

DATA
""" % lineage_cols)
		for s in self.samples:
			o.write("%s\t%s\n" % (s,lineage_cols.get(linresults[s]["main"],"#000000")))
		o.close()

		o = open(self.prefix+".dr.itol.txt","w")
		dr_cols = {"Sensitive":"#80FF00","Drug-resistant":"#00FFFF","MDR":"#8000FF","XDR":"#FF0000"}
		o.write("""DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL	Drug-Resistance
COLOR	#ff0000

LEGEND_TITLE	Drug resistance
LEGEND_SHAPES	1	1	1	1
LEGEND_COLORS	#80FF00	#00FFFF	#8000FF	#FF0000
LEGEND_LABELS	Sensitive	Drug-resisant	MDR	XDR

DATA
""")
		for s in self.samples:
			o.write("%s\t%s\n" % (s,dr_cols.get(results[s]["drtype"],"#000000")))
		o.close()
