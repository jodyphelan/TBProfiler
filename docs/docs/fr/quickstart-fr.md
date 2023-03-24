# Démarrage rapide

## Exemple de démarrage rapide

Exécuter l'ensemble du pipeline :

```
tb-profiler profile -1 /path/to/reads_1.fastq.gz -2 /path/to/reads_2.fastq.gz -p prefix
```

L'argument `-p` vous permet de fournir un préfixe aux fichiers de sortie résultants. Ceci est utile lorsque vous avez besoin d'exécuter plus d'un échantillon. Les fichiers BAM, VCF et resultat seront stockés dans les répertoires respectifs. Les résultats sont affichés au format json et texte.

## Exemple d'exécution

```
mkdir test_run; cd test_run
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR166/009/ERR1664619/ERR1664619_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR166/009/ERR1664619/ERR1664619_2.fastq.gz
tb-profiler profile -1 ERR1664619_1.fastq.gz -2 ERR1664619_2.fastq.gz -t 4 -p ERR1664619
cat results/ERR1664619.results.json
```

## Exécution à partir d'un fichier BAM existant

By using the `-a` option you can specify to use an existing BAM file instead of fastq files.

```
tb-profiler profile -a /path/to/bam -p test
```
! !! avertissement

    Les fichiers BAM ont dû être créés en utilisant la version du génome comme base de données qui peut être téléchargée ici. Ce génome a plusieurs numéros d'accès (ASM19595v2, NC_000962.3, GCF_000195955.2, etc...), ce qui peut prêter à confusion. Si vous pensez que votre référence est exactement la même séquence (la longueur devrait être de 4411532), vous pouvez créer une base de données avec le même nom de séquence que celui utilisé dans votre fichier BAM. Par exemple, si le nom de votre séquence est "NC_000962.3", vous pouvez le faire en exécutant ce qui suit :

    ```
    tb-profiler update_tbdb --match_ref /path/to/your/reference.fasta
    ```
