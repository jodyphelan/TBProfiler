# TBProfiler

This is the experimental commandline version of the TBProfiler described here: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-015-0164-0 and is available here: http://tbdr.lshtm.ac.uk

This repository contains a complete rewrite of the web version of TBProfiler. It allows the use of profiling through a command line interface and contains some additional functionality such as the ability to process minION data.

The pipeline aligns reads to the H37Rv reference BWA or minimap2 and then looks at the coverage across a number of different candidate regions. We also predict the number of reads supporting drug resistance variants as an insight into hetero-resistance (not applicable for minION data)

## Important changes as of v0.3.0
Version 0.3.0 features many changes to the code including the modularisation through creation of specific classes and functions for various functions. This will make it easier to maintain and add additional functionality. Support for minION has also been added. In the process some old functionality has not been added yet and some required arguments have changed. Please download v0.2.1 from the releases if you require the old code.

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

#### Quick start example
Run whole pipeline:
```
tb-profiler profile -1 /path/to/reads_1.fastq.gz -2 /path/to/reads_2.fastq.gz -p prefix
```
The prefix is usefull when you need to run more that one sample.
This will store BAM files, summary pileup and result files in respective directories.
Results are output in text format and json format.


Example run:
```
mkdir test_run; cd test_run
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR166/009/ERR1664619/ERR1664619_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR166/009/ERR1664619/ERR1664619_2.fastq.gz
../tb-profiler profile -1 ERR1664619_1.fastq.gz -2 ERR1664619_2.fastq.gz -t 4 -p ERR1664619
cat results/ERR1664619.results.txt
```

#### Running with an existing BAM file:

By using the ```-a``` option you can specify to use an existing BAM file instead of fastq files.
**Warning!!!**: The BAM files must have been created using the ensembl version of the genome which can be downloaded here:
```
ftp://ftp.ensemblgenomes.org/pub/release-32/bacteria//fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/dna/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.toplevel.fa.gz
```


The results from numerous runs can be collated into one table using the following command:
```
tb-profiler collate samples_file out_file
```
Where  ```samples_file``` is a list of prefixes of previously run samples and ```out_file``` is the name of the output file.


## Under the hood

The pipeline searches for small variants and big deletions associated with drug resistance. It will also report the lineage.

<img src="https://jodyphelan.github.io/img/TBProfiler.png">

### Adding new genes/mutations
To add new mutations navigate to the ```db``` directory and edit the ```drdb.txt``` file.
Add a new line corresponding to the desired variant with the following columns:

1. Drug - Drug name with no spaces
2. Genomic position - If more than one position affected, seperate with "/".
3. Reference nucleotides - String of nucleotides with length equal to the number of bases affected.
3. Alternate nucleotides - String of nucleotides with length equal to the number of bases affected.
4. Gene name
5. Mutation - String with the mutation

After editing the file run the ```parse_drdb.py <prefix>``` script using a prefix to generate a new database.
To use the database use the ```--db <prefix>``` option in ```tb-profiler```.

#### Examples:
##### A non-synonymous variant:

ETHIONAMIDE     1674484/1674485 AT      CC      inhA    Ile95Pro

##### A promoter mutation:

ISONIAZID       2156118 C       T       katG_promoter   C-7T

##### An indel:

PYRAZINAMIDE    2288953 CC      C       pncA    CC289C


## Change Log
#### v0.3
Modularisation of code into classes and functions
Support for minION

#### v0.2.1
* Collate data fix
* Generation of ITOL data files for visualisation

#### v0.2
* Allow for the choice of BWA or SNAP for Linux users
* Calling of deletions before small variant calling to avoid low quality variants around deletion breakpoints

#### v0.1
* Allow users to provide BAM file as input
* Ability to print out version


## Citation

If you would like to cite this work please use:

Coll, F. et al. Rapid determination of anti-tuberculosis drug resistance from whole-genome sequences. Genome Med. 5:51, (2015).
