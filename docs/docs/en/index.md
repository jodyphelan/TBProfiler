# TB-Profiler

TB-Profiler is a command-line tool used to process whole genome sequencing data to predict lineage and drug-resistance. The pipeline searches for small variants and big deletions associated with drug resistance. It will also report the lineage. By default it uses trimmomatic to trim the reads, BWA (or minimap2 for nanopore) to align to the reference genome and bcftools (or GATKv4/freebayes) to call variants. 

<img src="https://files.gitbook.com/v0/b/gitbook-legacy-files/o/assets%2F-M9cvGy4eVqvGN5UqFAr%2F-M9dZ5yUBJa3XVGD-wRl%2F-M9d_diTzl9Ae2kLM03R%2Ftb-profiler_uml.svg?alt=media&token=7cf9db30-a0c6-448b-a750-f4a041f34478">