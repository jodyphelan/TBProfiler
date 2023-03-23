# Execution sur plusieur echantillons

L'exécution du pipeline sur quelques échantillons est relativement simple, mais lorsque nous devons l'exécuter sur des centaines d'échantillons, il peut devenir un peu compliqué de taper chaque commande. Pour cela, il est possible d'utiliser la fonction `batch`.

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

Ici, vous pouvez fournir un fichier CSV avec les rubriques suivantes :

* id - ceci sera utilisé pour nommer les fichiers (obligatoire)
* read1 - la voie pour les reads forward
* read2 - a voie pour les reads reverse 
* bam - la voie pour les fichiers bam/cram 
* vcf - la voie pour les fichiers vcf 
* fasta - la voie pour les fichiers fasta

Chaque ligne doit contenir au moins le champ **id** et au moins l'un des champs du fichier d'entrée en fonction des données dont vous disposez.
