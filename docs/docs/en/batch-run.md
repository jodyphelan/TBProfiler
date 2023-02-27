# Run on many samples

To run the pipeline on a few samples if relatively straightforward, but when we need to run it on 100s of samples it can get a bit unwieldy to type every command out. For this purpose it is possible to use the `batch` function.

```
usage: tb-profiler batch [-h] --csv CSV [--args ARGS] [--jobs JOBS]
                         [--threads_per_job THREADS_PER_JOB] [--dir DIR]
                         [--temp TEMP] [--version]

optional arguments:
  -h, --help            show this help message and exit
  --csv CSV             CSV with samples and files (default: None)
  --args ARGS           Arguments to use with tb-profiler (default: None)
  --jobs JOBS, -j JOBS  Threads to use (default: 1)
  --threads_per_job THREADS_PER_JOB, -t THREADS_PER_JOB
                        Threads to use (default: 1)
  --dir DIR, -d DIR     Storage directory (default: .)
  --temp TEMP           Temp firectory to process all files (default: .)
  --version             show program's version number and exit
```

Here you can supply a CSV file with the following headings:

* id - this will be used to name the files (required)
* read1 - the path to the forward read
* read2 - the path to the reverse read 
* bam - the path to the bam/cram file 
* vcf - the path to the vcf file
* fasta - the path to the fasta file

Each line should have at least the **id** field and at least of of the input file fields depending on what data you have.