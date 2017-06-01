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

#### Quick start example
Run whole pipeline:
```
tb-profiler full -1 /path/to/reads_1.fastq.gz -2 /path/to/reads_2.fastq.gz -p prefix
```
The prefix is usefull when you need to run more that one sample.
This will store BAM files, summary pileup and result files in respective directories.
Results are output in text format and json format.

Example run:
```
mkdir test_run; cd test_run
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR166/009/ERR1664619/ERR1664619_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR166/009/ERR1664619/ERR1664619_2.fastq.gz
tb-profiler full -1 ERR1664619_1.fastq.gz -2 ERR1664619_2.fastq.gz -t 4 -p ERR1664619
cat results/ERR1664619.results.txt
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



## Citation

If you would like to cite this work please use:

Coll, F. et al. Rapid determination of anti-tuberculosis drug resistance from whole-genome sequences. Genome Med. 5:51, (2015).
