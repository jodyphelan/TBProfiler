import setuptools
from glob import glob

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
	data_files=[('share/tbprofiler',glob("db/*"))],
)
