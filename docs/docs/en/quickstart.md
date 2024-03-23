# Quickstart

## Quick start example

Run whole pipeline:

```
tb-profiler profile -1 /path/to/reads_1.fastq.gz -2 /path/to/reads_2.fastq.gz -p prefix
```

The `-p` argument allows you to provide a prefix to the resulting output files. This is useful when you need to run more that one sample. This will store BAM, VCF and result files in respective directories. Results are output in json and text format.

## Example run

```
mkdir test_run; cd test_run
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR166/009/ERR1664619/ERR1664619_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR166/009/ERR1664619/ERR1664619_2.fastq.gz
tb-profiler profile -1 ERR1664619_1.fastq.gz -2 ERR1664619_2.fastq.gz -t 4 -p ERR1664619
cat results/ERR1664619.results.json
```

## Running with an existing BAM file

By using the `-a` option you can specify to use an existing BAM file instead of fastq files. 

```
tb-profiler profile -a /path/to/bam -p test
```

!!! warning

    The BAM files must have been created using the version of the genome as the database which can be downloaded here. Confusingly, this genome has multiple accession numbers (ASM19595v2, NC_000962.3, GCF_000195955.2, etc...). If you believe your reference to be the exact same sequence (length should be 4411532) then you can create a database with the same sequence name as used in your BAM file. For example if your sequence name is "NC_000962.3" you can do this by executing the following:

    ```
    tb-profiler update_tbdb --match_ref /path/to/your/reference.fasta
    ```
