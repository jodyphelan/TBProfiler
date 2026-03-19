import os
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
	data_files=[
        (
            'share/tbprofiler/who_v2+',
            [x for x in glob("db/who_v2+/*") if not os.path.isdir(x)]
        ),
        (
            'share/tbprofiler/who_v2+/snpeff',
            [x for x in glob("db/who_v2+/snpeff/*") if not os.path.isdir(x)]
        ),
        (
            'share/tbprofiler/who_v2+/snpeff/data/Mycobacterium_tuberculosis_h37rv_tbprofiler',
            [x for x in glob("db/who_v2+/snpeff/data/Mycobacterium_tuberculosis_h37rv_tbprofiler/*") if not os.path.isdir(x)]
        ),
        (
            'share/tbprofiler/',
            [x for x in glob("db/*docx") if not os.path.isdir(x)]
        )
    ],
)
