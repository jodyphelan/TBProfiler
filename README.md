# TBProfiler

This is the experimental commandline version of the TBProfiler.
The pipeline aligns reads to the H37Rv reference using SNAP and then looks at the coverage across a number of different candidate regions.

#### Installation

```
git clone --recursive https://github.com/jodyphelan/TBProfiler.git
cd TBProfiler
bash install_prerequisites.sh 
echo "export PATH=\$PATH:$PWD" >> ~/.bashrc
```

#### Usage

Quick start
```
tb-profiler sv -1 /path/to/reads_1.fastq.gz -2/path/to/reads_2.fastq.gz -p prefix
```

