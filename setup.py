import setuptools


setuptools.setup(

	name="tbprofiler",

	version="3.0.5",
	packages=["tbprofiler","pathogenprofiler"],
	license="MIT",
	long_description="TBProfiler command line tool",
	scripts= [
		'tb-profiler'
		],
	data_files=[('share/tbprofiler',["db/tbdb.ann.txt","db/tbdb.barcode.bed","db/tbdb.bed","db/tbdb.dr.json","db/tbdb.fasta","db/tbdb.gff","db/tbdb.version.json","example_data/tbprofiler.test.fq.gz"])]
)
