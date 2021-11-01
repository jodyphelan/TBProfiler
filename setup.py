import setuptools

version = [l.strip() for l in open("tbprofiler/__init__.py") if "version" in l][0].split('"')[1]

setuptools.setup(

	name="tbprofiler",

	version=version,
	packages=["tbprofiler","pathogenprofiler"],
	license="MIT",
	long_description="TBProfiler command line tool",
	scripts= [
		'tb-profiler',
		'scripts/combine_vcf_variants.py',
		'scripts/rename_vcf_chrom.py',
		'scripts/add_dummy_AD.py'
		],
	data_files=[('share/tbprofiler',["db/tbdb.ann.txt","db/tbdb.barcode.bed","db/tbdb.bed","db/tbdb.dr.json","db/tbdb.fasta","db/tbdb.gff","db/tbdb.version.json","example_data/tbprofiler.test.fq.gz"])]
)
