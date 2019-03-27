#!/usr/bin/env sh
export CPP_INCLUDE_PATH=${CONDA_PREFIX}/include
export CXX_INCLUDE_PATH=${CONDA_PREFIX}/include
export CPLUS_INCLUDE_PATH=${CONDA_PREFIX}/include
export LIBRARY_PATH=${CONDA_PREFIX}/lib
conda install -y htslib libboost wget samtools bcftools bwa minimap2 bowtie2 tqdm parallel trimmomatic

cd /tmp/
wget https://github.com/dellytools/delly/archive/v0.8.1.tar.gz
tar -xvf v0.8.1.tar.gz
cd delly-0.8.1
make all
cp ./src/delly ${CONDA_PREFIX}/bin/delly_0.8.1
echo '#!/bin/sh' > ${CONDA_PREFIX}/bin/delly
echo export DYLD_FALLBACK_LIBRARY_PATH=${CONDA_PREFIX}/lib  >> ${CONDA_PREFIX}/bin/delly
echo ${CONDA_PREFIX}/bin/delly_0.8.1 '$@' >> ${CONDA_PREFIX}/bin/delly
chmod 755 ${CONDA_PREFIX}/bin/delly


cd /tmp/
wget https://github.com/jodyphelan/pathogen-profiler/archive/v1.1.tar.gz
tar -xvf v1.1.tar.gz
cd pathogen-profiler-1.1/
python setup.py install

cd /tmp/
wget https://github.com/jodyphelan/TBProfiler/archive/v2.1.tar.gz
tar -xvf v2.1.tar.gz
cd TBProfiler-2.1/
python setup.py install

cd /tmp/
mkdir ${CONDA_PREFIX}/share/tbprofiler
tb-profiler update_tbdb
rm -fr v0.8.1.tar.gz delly-0.8.1 tbdb TBProfiler-2.1 v2.1.tar.gz v1.1.tar.gz pathogen-profiler-1.1
