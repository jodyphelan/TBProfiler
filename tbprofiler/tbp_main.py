import json
import parse_csq
import db_compare
import mapping
import variant_calling
import files
import lineage
import deletions
import output

class tbp_seq_obj:
	params = {"fq1":False,"fq2":False,"bamfile":False,"mapping":True}
	def __init__(self,conf_file,prefix,fq1=False,fq2=False,bam=False,platform="Illumina",threads=1,outfmt="classic",db=False,stor_dir=".",verbose=False):
		if fq1 and files.filecheck(fq1) and files.verify_fq(fq1):
			self.params["fq1"] = fq1
		if fq2 and files.filecheck(fq2) and files.verify_fq(fq2):
			self.params["fq2"] = fq2
		if bam and files.filecheck(bam):
			self.params["bamfile"] = bam
			self.params["mapping"] = False
		else:
			self.params["bamfile"] = "%s/bam/%s.bam" % (stor_dir,prefix)
		self.params["stor_dir"] = stor_dir
		self.params["verbose"] = verbose
		self.params["outfmt"] = outfmt
		self.params["prefix"] = prefix
		self.params["dr_vcffile"] = "%s/vcf/%s.dr.vcf" % (stor_dir,prefix)
		self.params["vcffile"] = "%s/vcf/%s.vcf" % (stor_dir,prefix)
		self.params["gvcffile"] = "%s/vcf/%s.g.vcf.gz" % (stor_dir,prefix)
		self.params["temp_file"] = "%s.temp_file" % prefix
		self.params["temp_bam"] = "%s.temp.bam" % prefix
		self.params["temp_pileup"] = "%s.temp.pileup" % prefix
		self.params["platform"] = platform
		self.params["conf_file"] = conf_file
		self.params["depthfile"] = "%s.depth" % prefix
		self.params["txt_results"] = "%s/results/%s.results.txt" % (stor_dir,prefix)
		self.params["json_results"] = "%s/results/%s.results.json" % (stor_dir,prefix)
		self.params["threads"] = threads

		tmp = json.load(open(conf_file))
		for x in tmp:
			self.params[x] = tmp[x]
		if db:
			self.params["dr_json"] = "%s.json" % db
			self.params["dr_bed_file"] = "%s.bed" % db
			files.filecheck(self.params["dr_json"])
			files.filecheck(self.params["dr_bed_file"])
		if not fq1 and not fq1 and not bam:
			print "\nNo BAM or fastQ file provided. Please specify using -a or -1 (and -2).\n"
			quit()

	def call_genome_variants(self,gvcf=False):
		self.init_dirs()
		if self.params["mapping"]:
			self.map()
		else:
			files.index_bam(self)
		self.call_variants(gvcf=gvcf)

	def run_profiler(self):
		self.init_dirs()
		if self.params["mapping"]:
			self.map()
		else:
			files.index_bam(self)
		self.deletions()
		self.small_dr_variants()
		self.lineage()
		self.write_results()
		self.cleanup()

	def cleanup(self):
		if self.params["platform"]=="minION":
			files.rm_files(["%s%s" % (d,self.params["prefix"]) for d in [".temp.pileup",".temp.bam",".temp_file",".depth"]])
		else:
			files.rm_files(["%s%s" % (d,self.params["prefix"]) for d in [".temp.pileup",".temp.bam",".depth"]])

	def init_dirs(self):
		files.init_storage(self)

	def map(self):
		mapping.map(self)

	def small_dr_variants(self):
		variant_calling.call_variants(self,bed_file=self.params["dr_bed_file"],vcf_file=self.params["dr_vcffile"],caller="lofreq")
		variants = parse_csq.load_csq(self)
		self.small_dr_variants = db_compare.db_compare(self,variants)

	def call_variants(self,gvcf=False):
		variant_calling.call_variants(self,gvcf=gvcf)

	def lineage(self):
		self.lineage = lineage.lineage(self)

	def deletions(self):
		self.deletions = deletions.deletions(self)

	def write_results(self):
		if self.params["outfmt"] == "classic":
			output.write_results_2(self)
		elif self.params["outfmt"] == "new":
			output.write_results_1(self)
		else:
			print "Choose either 'new' or 'classic' as outfmt"
			quit()
