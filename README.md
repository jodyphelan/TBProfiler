# TBProfiler
[![Anaconda-Server Badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://bioconda.github.io/recipes/tb-profiler/README.html) [![Anaconda-Server Badge](https://img.shields.io/github/license/jodyphelan/TBProfiler.svg)](https://anaconda.org/bioconda/tb-profiler) [![Anaconda-Server Badge](https://img.shields.io/github/last-commit/jodyphelan/TBProfiler.svg)](https://github.com/jodyphelan/TBProfiler)

This repository contains a complete rewrite of the [web version of TBProfiler](http://tbdr.lshtm.ac.uk), described [here](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0650-x). It allows the use of profiling through a command line interface and contains some additional functionality such as the ability to process minION data.

The pipeline aligns reads to the H37Rv reference using bowtie2, BWA or minimap2 and then calls variants using SAMtools. These variants are then compared to a drug-resistance database. We also predict the number of reads supporting drug resistance variants as an insight into hetero-resistance (not applicable for minION data)

## Installation

##### Conda
Conda can function as a package manages are is available [here](https://docs.conda.io/en/latest/miniconda.html).
If you have conda make sure the bioconda  and conda-forge channels are added:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
Then you can install tb-profiler and all of its dependancies from the bioconda channel:
###### Linux
```
conda install -c bioconda tb-profiler
```
###### OSX
```
conda install -c bioconda tb-profiler samtools=1.9=h7c4ea83_11 ncurses=6.1=h0a44026_1002
```
##### Manually
It is possible to install manually. The following pre-requisites will be needed at runtime: *trimmomatic (>=v0.38), bwa (>=v0.7.17), minimap2 (>=v2.16), bowtie2 (>=v2.3.5), samtools (>=v1.9), bcftools (>=v1.9), GATK (>=v4.1.4.0), tqdm (>=v4.32.2) and parallel (>=v20190522)*. The pipeline should work and has been tested on the program versions indicated in parentheses.

To install both libraries run the following code:
```
pip3 install git+https://github.com/jodyphelan/TBProfiler.git
mkdir `python -c "import sys; print(getattr(sys, 'base_prefix', getattr(sys, 'real_prefix', sys.prefix)));"`
tb-profiler update_tbdb
```
You should then be able to run using ```tb-profiler```

## Usage

The first argument indicates the analysis type to perform. At the moment we currently only support the calling of small variants.

### Quick start example

Run whole pipeline:
```
tb-profiler profile -1 /path/to/reads_1.fastq.gz -2 /path/to/reads_2.fastq.gz -p prefix
```
The prefix is usefull when you need to run more that one sample. This will store BAM, VCF and result files in respective directories. Results are output in json and text format.

###### Example run
```
mkdir test_run; cd test_run
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR166/009/ERR1664619/ERR1664619_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR166/009/ERR1664619/ERR1664619_2.fastq.gz
tb-profiler profile -1 ERR1664619_1.fastq.gz -2 ERR1664619_2.fastq.gz -t 4 -p ERR1664619
cat results/ERR1664619.results.json
```

###### Running with an existing BAM file

By using the -a option you can specify to use an existing BAM file instead of fastq files. **Warning!!!**: The BAM files must have been created using the version of the genome as the database which can be downloaded [here](ftp://ftp.ensemblgenomes.org/pub/release-32/bacteria//fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/dna/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.toplevel.fa.gz). Confusingly, this genome has multiple accession numbers (ASM19595v2,NC_000962.3,GCF_000195955.2, etc...). If you believe your reference to be the exact same sequence (length should be 4411532) then you can create a database with the same sequence name as used in your BAM file. For example if your sequence name is "NC_000962.3" you can do this by either:
1. Creating the new database files using the `--seqname NC_000962.3` option from the `parse_db.py` script in the [tbdb repo](https://github.com/jodyphelan/tbdb). Then loading it using `tb-profiler load_library /path/to/lib`.
2. Or applying a quick fix to replace references to "Chromosome" in all existing database files e.g using ```sed -i 's/Chromosome/NC_000962.3/' `python -c "import sys;print(getattr(sys, 'base_prefix', getattr(sys, 'real_prefix', sys.prefix)))"`/share/tbprofiler/tbdb* && samtools faidx `python -c "import sys;print(getattr(sys, 'base_prefix', getattr(sys, 'real_prefix', sys.prefix)))"`/share/tbprofiler/tbdb.fasta```

###### Summarising runs

The results from numerous runs can be collated into one table using the following command:
```
tb-profiler collate
```
This will automatically create a number of colled result files from all the individual result files in the *result* directory. If you would like to generate this file for a subset of the runs you can provide a list with the run sames using the `--samples` flag. The prefix for the output files is *tbprofiler* by default but this can be changed with the `--prefix` flag.

## Mutation database
TBProfiler ships with a default database. The development of the mutation library is hosted on the [tbdb repository](https://github.com/jodyphelan/tbdb). Please visity this repo if you would like to get involved in the database or would like to modify and create your own.

If you would like to use an altered database you can load the config file produced by `parse_db.py` as such:
```
tb-profiler load_library [config.json]
```

### Non-H37Rv databases
It is possible run TBProfiler on another reference genome. Although there is currently no helper tool to create the databases for other references automatically, checkout the [tbdb repository](https://github.com/jodyphelan/tbdb) to find out more about what you need.

## Under the hood
The pipeline searches for small variants and big deletions associated with drug resistance. It will also report the lineage.
By default it uses Trimmomatic to trim the reads, BWA (or minimap2 for nanopore) to align to the reference genome and GATK (open source v4) to call variants.

<img src="https://github.com/jodyphelan/jodyphelan.github.io/raw/master/TBProfiler_pipeline.png">

## ITOL files
Several files are produced by the `tb-profile collate` function. Among these are several config files that can be used with [iTOL](http://itol.embl.de/) to annotate phylogenetic trees. A small tree and config files have been placed in the *example_data* directory. To use navigate to the iTOL website and upload the *tbprofiler.tree* file using the upload button on the navigation bar. Once this has been uploaded you will be taken to a visualisation of the tree. To add the annotation, click on the '+' button on the lower right hand corner and select the iTOL config files. You should now see a figure similar to the one below. The following annotations are included:
 - Lineage
 - Drug resistance classes (Sensitive, drug-resistant, MDR, XDR)
 - Drug resistance calls for individual drugs, were filled circles represent resistance.

<img src="https://github.com/jodyphelan/jodyphelan.github.io/raw/master/img/itol_example.png">

## Issues
Please raise them using the [Issues](https://github.com/jodyphelan/TBProfiler/issues) page.

## FAQ

Will populate this once we get some frequently asked questions!

## Future plans
 - Add in capability to perform basic phylogenetic functions
 - Add in levels of resistance to mutations
