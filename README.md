![TB Profiler Logo](https://raw.githubusercontent.com/jodyphelan/jodyphelan_github_io_old/master/img/tb-profiler-logo-rectangle.png)

[![Anaconda-Server Badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://bioconda.github.io/recipes/tb-profiler/README.html) [![Anaconda-Server Badge](https://img.shields.io/github/license/jodyphelan/TBProfiler.svg)](https://anaconda.org/bioconda/tb-profiler) [![Anaconda-Server Badge](https://img.shields.io/github/last-commit/jodyphelan/TBProfiler.svg)](https://github.com/jodyphelan/TBProfiler) [![Codacy Badge](https://app.codacy.com/project/badge/Grade/1f1e0523018c4871a3b56583f59beae9)](https://www.codacy.com/gh/jodyphelan/TBProfiler/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=jodyphelan/TBProfiler&amp;utm_campaign=Badge_Grade) [![Codacy Badge](https://app.codacy.com/project/badge/Coverage/1f1e0523018c4871a3b56583f59beae9)](https://www.codacy.com/gh/jodyphelan/TBProfiler/dashboard?utm_source=github.com&utm_medium=referral&utm_content=jodyphelan/TBProfiler&utm_campaign=Badge_Coverage)



This repository contains a complete rewrite of the [web version of TBProfiler](http://tbdr.lshtm.ac.uk), described [here](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0650-x). It allows the use of profiling through a command line interface and contains some additional functionality such as the ability to process minION data.

The pipeline aligns reads to the H37Rv reference using bowtie2, BWA or minimap2 and then calls variants using bcftools. These variants are then compared to a drug-resistance database. We also predict the number of reads supporting drug resistance variants as an insight into hetero-resistance \(not applicable for minION data\)

## Documentation

This page has all the info you need to get started, however additional (and more organised!) documentation is available [here](https://jodyphelan.gitbook.io/tb-profiler/). We have also have some basic translation to other languages:  [:brazil:](https://jodyphelan.gitbook.io/tb-profiler/translations/portugues)[:netherlands:](https://jodyphelan.gitbook.io/tb-profiler/translations/nederlands). Please contact us if you would like to improve a translation or add a new one!

## Keeping up to date
TBProfiler is under constant rapid development. If you plan to use the program in your work please make sure you are using the most up to date version! Similarly, the database is not static and is continuously being improved so make sure you are using the most latest version. If you use TBProfiler in your work please state the version of both the tool and the database as they are deveoped independantly from each other.

## Installation

### Conda

Conda can function as a package manages are is available [here](https://docs.conda.io/en/latest/miniconda.html). If you have conda make sure the bioconda and conda-forge channels are added:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Then you can install tb-profiler and all of its dependancies from the bioconda channel:

#### Linux

```bash
conda install -c bioconda tb-profiler
```

#### OSX

```bash
conda install -c bioconda tb-profiler samtools=1.9=h7c4ea83_11 ncurses=6.1=h0a44026_1002
```

### Manually

It is possible to install manually. The following pre-requisites will be needed at runtime: _trimmomatic \(&gt;=v0.38\), bwa \(&gt;=v0.7.17\), minimap2 \(&gt;=v2.16\), bowtie2 \(&gt;=v2.3.5\), samtools \(&gt;=v1.10\), bcftools \(&gt;=v1.10\), GATK \(&gt;=v4.1.4.0\), freebayes \(&gt;=v1.3.2\), tqdm \(&gt;=v4.32.2\) parallel \(&gt;=v20190522\) samclip \(&gt;=v0.4.0\)_ and snpEff \(&gt;=v5.0.0\)_. The pipeline should work and has been tested on the program versions indicated in parentheses.

To install the library run the following code:
```bash
pip3 install git+https://github.com/jodyphelan/TBProfiler.git
pip3 install git+https://github.com/jodyphelan/pathogen-profiler.git
mkdir `python -c "import sys; print(getattr(sys, 'base_prefix', getattr(sys, 'real_prefix', sys.prefix)));"`
tb-profiler update_tbdb
```

You should then be able to run using `tb-profiler`

## Usage

The first argument indicates the analysis type to perform. At the moment we currently only support the calling of small variants.

### Quick start example

Run whole pipeline:

```text
tb-profiler profile -1 /path/to/reads_1.fastq.gz -2 /path/to/reads_2.fastq.gz -p prefix
```

The prefix is usefull when you need to run more that one sample. This will store BAM, VCF and result files in respective directories. Results are output in json and text format.

#### Example run

```bash
mkdir test_run
cd test_run
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR166/009/ERR1664619/ERR1664619_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR166/009/ERR1664619/ERR1664619_2.fastq.gz
tb-profiler profile -1 ERR1664619_1.fastq.gz -2 ERR1664619_2.fastq.gz -t 4 -p ERR1664619 --txt
cat results/ERR1664619.results.txt
```

#### Spoligotyping

Experimental spoligotyping can be performed by adding `--spoligotype` to the command. This is enabled for bam and fastq input.

#### Running with an existing BAM file

By using the -a option you can specify to use an existing BAM file instead of fastq files. **Warning!!!**: The BAM files must have been created using the version of the genome as the database which can be downloaded [here](ftp://ftp.ensemblgenomes.org/pub/release-32/bacteria//fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/dna/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.toplevel.fa.gz). Confusingly, this genome has multiple accession numbers \(ASM19595v2,NC\_000962.3,GCF\_000195955.2, etc...\). If you believe your reference to be the exact same sequence \(length should be 4411532\) then you can create a database with the same sequence name as used in your BAM file. For example if your sequence name is "NC\_000962.3" you can run

```bash
tb-profiler update_tbdb --seqname NC_000962.3
```

### Summarising runs

The results from numerous runs can be collated into one table using the following command:

```bash
tb-profiler collate
```

This will automatically create a number of colled result files from all the individual result files in the _result_ directory. If you would like to generate this file for a subset of the runs you can provide a list with the run sames using the `--samples` flag. The prefix for the output files is _tbprofiler_ by default but this can be changed with the `--prefix` flag.

### Writing your own summary scripts

The `collate` function extracts the drug-resistance mutations and lineage, however you may want to extract more features that are present in the individual json result files. I have created a little tutorial on how to do this [here](https://jodyphelan.gitbook.io/tb-profiler/writing-a-custom-collate-script).

## Mutation database

TBProfiler ships with a default database. The development of the mutation library is hosted on the [tbdb repository](https://github.com/jodyphelan/tbdb). Please visity this repo if you would like to get involved in the database or would like to modify and create your own.

If you would like to use an altered database you can download the tbdb repo, make the required changes and run the following code from within the tbdb repo directory:

```text
tb-profiler create_db --prefix <new_library_name>
tb-profiler load_library --prefix <new_library_name>
```

### Non-H37Rv databases

It is possible run TBProfiler on another reference genome. Although there is currently no helper tool to create the databases for other references automatically, checkout the [tbdb repository](https://github.com/jodyphelan/tbdb) to find out more about what you need.

## Under the hood

The pipeline searches for small variants and big deletions associated with drug resistance. It will also report the lineage. By default it uses Trimmomatic to trim the reads, BWA \(or minimap2 for nanopore\) to align to the reference genome and bcftools (or GATKv4/freebayes) to call variants. Here is an example schematic for illumina paried end data:

![](https://github.com/jodyphelan/jodyphelan.github.io/raw/master/tb-profiler_uml.svg)

## ITOL files

Several files are produced by the `tb-profile collate` function. Among these are several config files that can be used with [iTOL](http://itol.embl.de/) to annotate phylogenetic trees. A small tree and config files have been placed in the _example\_data_ directory. To use navigate to the iTOL website and upload the _tbprofiler.tree_ file using the upload button on the navigation bar. Once this has been uploaded you will be taken to a visualisation of the tree. To add the annotation, click on the '+' button on the lower right hand corner and select the iTOL config files. You should now see a figure similar to the one below. The following annotations are included:

*   Lineage
*   Drug resistance classes \(Sensitive, drug-resistant, MDR, XDR\)
*   Drug resistance calls for individual drugs, were filled circles represent resistance.

![](https://github.com/jodyphelan/jodyphelan.github.io/raw/master/img/itol_example.png)

## Citation

If you plan to use TB-Profiler in your work please cite [the paper](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0650-x) and include both the version of `tb-profiler` and the database.

## Issues

Please raise them using the [Issues](https://github.com/jodyphelan/TBProfiler/issues) page. Please see the [contributing guidelines](https://github.com/jodyphelan/TBProfiler/blob/master/docs/CONTRIBUTING.md) and the [code of conduct](https://github.com/jodyphelan/TBProfiler/blob/master/CODE_OF_CONDUCT.md).

## FAQ
Q. How can I use tb-profiler to produce a phylogenetic tree?

A. At the moment TB-Profiler does not provide any functionality to actually create the tree. The config files for itol are generated, however the tree must be created with other software. It is planned for future releases, however it might take some time to be fully developed. I would suggest alternative software for the phylogenetic reconstruciton step (such as snippy + iqtree)

## Future plans

Please see the [roadmap](https://github.com/jodyphelan/TBProfiler/projects/1) to view future plans as well as features currently in development
