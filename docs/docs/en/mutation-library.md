# Mutation library

TB-Profiler ships with a default database. The development of the mutation library is hosted on the [tbdb repository](https://github.com/jodyphelan/tbdb). Please visit this repo if you would like to get involved in the database or would like to modify and create your own.

## Why is there a separate github repository?

With analysis pipelines pretty much standardised, it is evident that accuracy of prediction is affected mostly by the underlying library of mutations. As new evidence for the inclusion or exclusion of mutations is generated, there is a constant need to update and re-evaluate the mutation library. Moreover, it is important for the control of the library to be put in the hands of the end-users. By hosting the library on a separate repository (rather than buried deep in the profiling tool code) it makes it easier to find out exactly which mutations are present. Additionally, github has a number of useful features which can be utilised:

* All changes to the library are tracked and can be easily be reverted to previous versions.
* Users can raise concerns or discuss the library using the [Issues](https://github.com/jodyphelan/tbdb/issues) section of the github repository.
* Different versions of the library can be maintained using [Forks](https://help.github.com/en/articles/fork-a-repo). Allowing users to experiment with the library without affecting the main project. These changes can then be merged into the main repo after the changes are reviewed.
* Multiple users/developers can contribute towards the library.

!!! TLDR 
    In a nutshell: hosting it separately makes it easier to update the library.

# Want to contribute?
f you think a mutation should be removed or added please raise and issue [here](https://github.com/jodyphelan/tbdb/issues). If you want to help curate the library, leave a comment [here](https://github.com/jodyphelan/tbdb/issues/4).

## Adding/removing mutations
Mutations can be added by submitting a pull request on a branch modified tbdb.csv file. If that previous sentence made no sense to you then you can suggest a change using an [issue](https://github.com/jodyphelan/tbdb/issues) and we will try help.

## How does it work?

The mutations are listed in [tbdb.csv](https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv). These are parsed by `tb-profiler create_db` to generate the json formatted database used by TB-Profiler along with a few more files. Mutations can be removed and added from tbdb.csv and a new library can be built using `tb-profiler create_db`.

## tbdb.csv
This is a CSV file which must contain the following column headings:
1. Gene - These can be the gene names (e.g. rpoB) or locus tag (e.g. Rv0667).
2. Mutation - These must follow the [hgvs nomenclature](http://varnomen.hgvs.org/). More info on this down below.
3. Drug - Name of the drug
4. Confers - This should be set to 'resistance'
5. Interaction - This should be set to empty
6. Literature - Any literature which provides evidence for the mutation. I would recommend using pubmed IDs (e.g. PMC3315572) but in theory anything can be put here. Multiple entries can be separated with ";".
7. WHO Confidence - The confidence as given by the WHO mutation catalogue

The first three columns must contain a value, however literature may remain empty. Additional columns may be added and will be built into the json library, and can be output in the tb-profiler results.

## Mutation format
Mutations must follow the HGVS nomenclature. Information on this format can be found [here](http://varnomen.hgvs.org/). The following types of mutations are currently allowed:

* Amino acid substitutions. Example: S450L in rpoB would be p.Ser450Leu
* Deletions in genes. Example: Deletion of nucleotide 758 in tlyA would be c.758del
* Insertion in genes. Example: Insertion of GT between nucleotide 1850 and 1851 in katG would be c.1850_1851insGT
* SNPs in non-coding RNAs. Example: A to G at position 1401 in rrs would be n.1401a>g
* SNPs in gene promoters. Example: A to G 7 bases 5' of the start codon in pncA c.-7A>G

!!! warning

    The mutations and resulting library files are in reference to the H37Rv (NC_000962.3/AL123456.3) reference genome**

## Confidence values

The confidence values are taken from the WHO catalogue but this can be changed to whatever you want.

# Generating a new library

Just download the repository using `git clone https://github.com/jodyphelan/tbdb.git`. This will generate a folder with all the required files. Then you can run `tb-profiler create_db`

```
usage: tb-profiler create_db [-h] [--prefix PREFIX] [--csv CSV] [--watchlist WATCHLIST] [--spoligotypes SPOLIGOTYPES]
                             [--spoligotype_annotations SPOLIGOTYPE_ANNOTATIONS] [--barcode BARCODE] [--bedmask BEDMASK]
                             [--amplicon_primers AMPLICON_PRIMERS] [--match_ref MATCH_REF]
                             [--other_annotations OTHER_ANNOTATIONS] [--custom] [--db_name DB_NAME] [--db_commit DB_COMMIT]
                             [--db_author DB_AUTHOR] [--db_date DB_DATE] [--include_original_mutation] [--load]
                             [--no_overwrite] [--dir DIR] [--temp TEMP] [--version]

optional arguments:
  -h, --help            show this help message and exit
  --prefix PREFIX, -p PREFIX
                        The input CSV file containing the mutations (default: tbdb)
  --csv CSV, -c CSV     The prefix for all output files (default: tbdb.csv)
  --watchlist WATCHLIST, -w WATCHLIST
                        A csv file containing genes to profile but without any specific associated mutations (default:
                        tbdb.watchlist.csv)
  --spoligotypes SPOLIGOTYPES
                        A file containing a list of spoligotype spacers (default: spoligotype_spacers.txt)
  --spoligotype_annotations SPOLIGOTYPE_ANNOTATIONS
  --barcode BARCODE     A bed file containing lineage barcode SNPs (default: barcode.bed)
  --bedmask BEDMASK     A bed file containing a list of low-complexity regions (default: mask.bed)
  --amplicon_primers AMPLICON_PRIMERS
                        A file containing a list of amplicon primers (default: None)
  --match_ref MATCH_REF
                        The prefix for all output files (default: None)
  --other_annotations OTHER_ANNOTATIONS
                        A CSV containing gene, mutation, drug and confidence columns (default: tbdb.other_annotations.csv)
  --custom              Tells the script this is a custom database, this is used to alter the generation of the version file
                        (default: False)
  --db_name DB_NAME     Overrides the name of the database in the version file (default: None)
  --db_commit DB_COMMIT
                        Overrides the commit string of the database in the version file (default: None)
  --db_author DB_AUTHOR
                        Overrides the author of the database in the version file (default: None)
  --db_date DB_DATE     Overrides the date of the database in the version file (default: None)
  --include_original_mutation
                        Include the original mutation (before reformatting) as part of the variant annotaion (default: False)
  --load                Automaticaly load database (default: False)
  --no_overwrite        Don't load if existing database with prefix exists (default: False)
  --dir DIR, -d DIR     Storage directory (default: .)
  --temp TEMP           Temp firectory to process all files (default: .)
  --version             show program's version number and exit
  ```

If you run it without any arguments it will generate the database files with the prefix `tbdb`. This can be changed by using the `--prefix` argument.

## Using alternate reference names

The tbdb database will assume you have mapped your data to a reference with "Chromosome" as the sequence name. If your reference sequence is the same but has a differenct name e.g "NC_000962.3". You can generate a database with an alternate sequence name using the `--match_ref /path/to/your/reference.fasta` flag.

## Watchlist

There are some genes it may be of interest to record mutations even if we do not have any specific associated mutaitons. To allow this funcitonality we have included a "watchlist" file. To include genes just add them and the associated drug(s) to the tbdb.watchlist.csv file.
