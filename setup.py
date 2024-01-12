import setuptools

version = [l.strip() for l in open("tbprofiler/__init__.py") if "version" in l][0].split('"')[1]

setuptools.setup(

	name="tbprofiler",

	version=version,
	packages=["tbprofiler"],
	license="GPLv3",
	long_description="TBProfiler command line tool",
	scripts= [
		'tb-profiler',
		'scripts/tb-profiler-tools'
		],
	data_files=[('share/tbprofiler',["db/default_template.docx","db/tbdb.mask.bed","db/tbdb.spoligotype_spacers.txt","db/tbdb.spoligotype_list.csv","db/tbdb.barcode.bed","db/tbdb.bed","db/tbdb.dr.json","db/tbdb.fasta","db/tbdb.gff","db/tbdb.version.json","db/tbdb.variables.json","example_data/tbprofiler.test.fq.gz"])]
)
