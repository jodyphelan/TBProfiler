# TBProfiler

This is the experimental commandline version of the TBProfiler.
The pipeline aligns reads to the H37Rv reference using SNAP and then looks at the coverage across a number of different candidate regions.



### Installation

```
git clone --recursive https://github.com/jodyphelan/TBProfiler.git
cd TBProfiler
bash install_prerequisites.sh 
echo "export PATH=\$PATH:$PWD" >> ~/.bashrc
```

For OSX use:
```
bash osX_install_prerequisites.sh
```
### Usage

The first argument indicates the analysis type to perform. At the moment we currently only support the calling of small variants using ```sv``` or the detection of large deletions using ```del```.
#### Quick start 
Run whole pipeline:
```
tb-profiler full -1 /path/to/reads_1.fastq.gz -2 /path/to/reads_2.fastq.gz -p prefix
```

#### Step by Step
Map reads and perform pileup:
```
tb-profiler mapping -1 /path/to/reads_1.fastq.gz -2 /path/to/reads_2.fastq.gz -p prefix
``` 
Look for small variations (SNPs and small INDELs) potentially causing drug resistance:
```
tb-profiler sv -p prefix
```
Look for large deletions 
```
tb-profiler del -p prefix
```

