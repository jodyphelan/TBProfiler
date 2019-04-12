import pathogenprofiler as pp
import json
def profile_vcf(filename,conf):
	params = conf.copy()
	params["tmpvcf"] = pp.get_random_file(extension=".vcf.gz")
	params["tmpcsq"] = pp.get_random_file(extension=".vcf.gz")
	params["filename"] = filename
	params["tmphdr"] = pp.get_random_file()
	params["tmptxt"] = pp.get_random_file()
	l=""
	for l in pp.cmd_out("bcftools view %(filename)s -h | grep \"^##FORMAT=<ID=AD\"" % params):
		pass
	AD_found = False if l=="" else True
	if AD_found==False:
		open(params["tmphdr"],"w").write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n")
		pp.run_cmd("bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t.[\\t0,100]\\n' %(filename)s > %(tmptxt)s" % params)
		pp.run_cmd("bgzip %(tmptxt)s" % params)
		pp.run_cmd("tabix -s 1 -b 2 -p vcf %(tmptxt)s.gz" % params)
		pp.run_cmd("bcftools view -a %(filename)s | bcftools annotate -a %(tmptxt)s.gz -c CHROM,POS,REF,ALT,-,FMT/AD -h %(tmphdr)s -Oz -o %(tmpvcf)s" % params)
	else:
		pp.run_cmd("bcftools view -a %(filename)s -Oz -o %(tmpvcf)s" % params)
	pp.run_cmd("bcftools view -T %(bed)s %(tmpvcf)s | bcftools csq -f %(ref)s -g %(gff)s  -Oz -o %(tmpcsq)s -p a" % params)
	csq_bcf_obj = pp.bcf(params["tmpcsq"])
	csq = csq_bcf_obj.load_csq(ann_file=conf["ann"])
	results = {"variants":[],"missing_pos":[],"qc":{"pct_reads_mapped":"NA","num_reads_mapped":"NA"}}
	for sample in csq:
		results["variants"]  = csq[sample]
	all_bcf_obj = pp.bcf(params["tmpvcf"])
	mutations = all_bcf_obj.get_bed_gt(conf["barcode"],conf["ref"])
	if "C" in mutations["Chromosome"][325505] and  mutations["Chromosome"][325505]["C"]==50:  mutations["Chromosome"][325505] = {"T":25}
	if "G" in mutations["Chromosome"][599868] and  mutations["Chromosome"][599868]["G"]==50:  mutations["Chromosome"][599868] = {"A":25}
	if "C" in mutations["Chromosome"][931123] and  mutations["Chromosome"][931123]["C"]==50:  mutations["Chromosome"][931123] = {"T":25}
	if "T" in mutations["Chromosome"][1759252] and  mutations["Chromosome"][1759252]["T"]==50:  mutations["Chromosome"][1759252] = {"G":25}
	json.dump(mutations,open("dump.json","w"))
	barcode_mutations = pp.barcode(mutations,conf["barcode"])
	results["barcode"] = barcode_mutations
	results = pp.db_compare(db_file=conf["json_db"],mutations=results)
	bed_regions = pp.load_bed(conf["bed"],[4],4)
	missing_regions = {gene:"NA" for gene in bed_regions}
	results["missing_regions"] = missing_regions
	if AD_found:
		pp.run_cmd("rm %(tmpcsq)s" % params)
	else:
		pp.run_cmd("rm %(tmpcsq)s %(tmphdr)s %(tmptxt)s*" % params)
	return results
