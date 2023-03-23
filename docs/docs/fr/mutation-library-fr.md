# Bibliothèque de mutation

TB-Profiler est livré avec une base de données par défaut. Le développement de la bibliothèque de mutation est hébergé sur le [tbdb repository](https://github.com/jodyphelan/tbdb). Veuillez visiter ce dossier si vous souhaitez vous impliquer dans la base de données ou si vous souhaitez modifier et créer votre propre base de données.

## Pourquoi y a-t-il un repositoire github séparé ?

Les pipelines d'analyse étant à peu près standardisés, il est évident que la précision de la prédiction est principalement affectée par la bibliothèque de mutations sous-jacente. Au fur et à mesure que de nouvelles preuves d'inclusion ou d'exclusion de mutations sont produites, il est constamment nécessaire de mettre à jour et de réévaluer la bibliothèque de mutations. En outre, il est important que le contrôle de la bibliothèque soit confié aux utilisateurs finaux. En hébergeant la bibliothèque sur un repositoire séparé (plutôt qu'enfouie dans le code de l'outil de profilage), il est plus facile de savoir exactement quelles mutations sont présentes. En outre, github dispose d'un certain nombre de fonctionnalités utiles qui peuvent être utilisées :

* Toutes les modifications apportées à la bibliothèque sont suivies et peuvent être facilement ramenées à des versions antérieures.
* Les utilisateurs peuvent soulever des problèmes ou discuter de la bibliothèque en utilisant la section [Issues](https://github.com/jodyphelan/tbdb/issues) du repositoire github.
* Différentes versions de la bibliothèque peuvent être maintenues en utilisant [Forks](https://help.github.com/en/articles/fork-a-repo). Cela permet aux utilisateurs d'expérimenter avec la bibliothèque sans affecter le projet principal. Ces changements peuvent ensuite être fusionnés dans le repo principal après avoir été revus.
* Plusieurs utilisateurs/développeurs peuvent contribuer à la bibliothèque.

! !! TLDR 
    En bref : l'hébergement séparé facilite la mise à jour de la bibliothèque.

# Vous voulez contribuer ?
Si vous pensez qu'une mutation devrait être supprimée ou ajoutée, soulevez un problème [ici](https://github.com/jodyphelan/tbdb/issues). Si vous souhaitez contribuer à l'enrichissement de la bibliothèque, laissez un commentaire [ici](https://github.com/jodyphelan/tbdb/issues/4).

## Ajouter/supprimer des mutations
Les mutations peuvent être ajoutées en soumettant une pull request sur une branche modifiée du fichier tbdb.csv. Si la phrase précédente n'a aucun sens pour vous, vous pouvez suggérer un changement en utilisant un [issue](https://github.com/jodyphelan/tbdb/issues) et nous essayerons de vous aider.

## Comment cela fonctionne-t-il ?

Les mutations sont listées dans [tbdb.csv](https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv). Elles sont analysées par `tb-profiler create_db` pour générer la base de données au format json utilisée par TB-Profiler ainsi que quelques autres fichiers. Les mutations peuvent être supprimées et ajoutées à partir de tbdb.csv et une nouvelle bibliothèque peut être construite en utilisant `tb-profiler create_db`.

## tbdb.csv
Il s'agit d'un fichier CSV qui doit contenir les titres de colonnes suivants :
1. Gène - Il peut s'agir du nom du gène (par exemple, rpoB) ou de l'étiquette du locus (par exemple, Rv0667).
2. Mutation - Elle doit être conforme à la [nomenclature hgvs](http://varnomen.hgvs.org/). Plus d'informations à ce sujet ci-dessous.
3. Drug (médicament) - Nom du médicament
4. Confers - Ce paramètre doit être réglé sur "résistance
5. Interaction - Ce champ doit être vide
6. Littérature - Toute littérature fournissant des preuves de la mutation. Je recommande d'utiliser les ID de Pubmed (par exemple PMC3315572), mais en théorie, tout peut être mis ici. Les entrées multiples peuvent être séparées par des " ;".
7. WHO Confidence - La confiance donnée par le catalogue des mutations de l'OMS.

Les trois premières colonnes doivent contenir une valeur, mais la littérature peut rester vide. Des colonnes supplémentaires peuvent être ajoutées et seront intégrées dans la bibliothèque json, et peuvent être affichées dans les résultats de tb-profiler.

## Format de la mutation
Les mutations doivent suivre la nomenclature HGVS. Des informations sur ce format sont disponibles [ici](http://varnomen.hgvs.org/). Les types de mutations suivants sont actuellement autorisés :

* Substitutions d'acides aminés. Exemple : S450L dans rpoB serait p.Ser450Leu
* Délétions dans les gènes. Exemple : La délétion du nucléotide 758 dans tlyA serait c.758del.
* Insertion dans les gènes. Exemple : L'insertion de GT entre les nucléotides 1850 et 1851 dans le gène katG serait c.1850_1851insGT.
* SNP dans les ARN non codants. Exemple : L'insertion de A à G à la position 1401 dans l'ARS serait n.1401a>g.
* SNP dans les promoteurs de gènes. Exemple : A à G à 7 bases 5' du codon de départ dans pncA c.-7A>G

! !! avertissement

    Les mutations et les fichiers de bibliothèques qui en résultent se réfèrent au génome de référence H37Rv (NC_000962.3/AL123456.3)**

## Valeurs de confiance

Les valeurs de confiance sont tirées du catalogue de l'OMS, mais elles peuvent être modifiées à votre guise.

# Génération d'une nouvelle bibliothèque

Téléchargez simplement le repositoire en utilisant `git clone https://github.com/jodyphelan/tbdb.git`. Cela va générer un dossier avec tous les fichiers nécessaires. Ensuite, vous pouvez lancer `tb-profiler create_db`

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

Si vous l'executez sans aucun argument, il générera les fichiers de base de données avec le préfixe `tbdb`. Ceci peut être modifié en utilisant l'argument `--prefix`.

## Utilisation de noms de référence alternatifs

La base de données tbdb suppose que vous avez mis en correspondance vos données avec une référence dont le nom de séquence est "Chromosome". Si votre séquence de référence est la même mais a un nom différent, par exemple "NC_000962.3". Vous pouvez générer une base de données avec un nom de séquence alternatif en utilisant l'option `--match_ref /path/to/your/reference.fasta`.

## Liste de surveillance

Pour certains gènes, il peut être intéressant d'enregistrer les mutations même si nous n'avons pas de mutations spécifiques associées. Pour permettre cette fonctionnalité, nous avons inclus un fichier "watchlist". Pour inclure des gènes, il suffit de les ajouter, ainsi que le(s) médicament(s) associé(s), au fichier tbdb.watchlist.csv.
