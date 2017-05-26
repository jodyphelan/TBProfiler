# TBProfiler

This is the experimental commandline version of the TBProfiler described here: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-015-0164-0 and is available here: http://tbdr.lshtm.ac.uk

The pipeline aligns reads to the H37Rv reference using SNAP and then looks at the coverage across a number of different candidate regions.

## Installation

```
git clone --recursive https://github.com/jodyphelan/TBProfiler.git
cd TBProfiler
bash install_prerequisites.sh
echo "export PATH=$PWD:$PATH" >> ~/.bashrc
```

For OSX use:
```
bash osX_install_prerequisites.sh
```
## Usage

The first argument indicates the analysis type to perform. At the moment we currently only support the calling of small variants or the detection of large deletions.

#### Quick start
Run whole pipeline:
```
tb-profiler full -1 /path/to/reads_1.fastq.gz -2 /path/to/reads_2.fastq.gz -p prefix
```
The prefix is usefull when you need to run more that one sample.
This will store BAM files, summary pileup and result files in respective directories.
Results are output in text format and json format.

The results from numerous runs can be collated into one table using the following command:
```
tb-profiler collate samples_file out_file
```
Where  ```samples_file``` is a list of prefixes of previously run samples and ```out_file``` is the name of the output file.


## Under the hood

The pipeline searches for small variants and big deletions associated with drug resistance. It will also report the lineage.

### Adding new genes/mutations
This feature is coming soon

## Citation

If you would like to cite this work please use:

Coll, F. et al. Rapid determination of anti-tuberculosis drug resistance from whole-genome sequences. Genome Med. 5:51, (2015).
