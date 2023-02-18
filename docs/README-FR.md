![TB Profiler Logo](https://raw.githubusercontent.com/jodyphelan/jodyphelan_github_io_old/master/img/tb-profiler-logo-rectangle.png)

[![Anaconda-Server Badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://bioconda.github.io/recipes/tb-profiler/README.html) [![Anaconda-Server Badge](https://img.shields.io/github/license/jodyphelan/TBProfiler.svg)](https://anaconda.org/bioconda/tb-profiler) [![Anaconda-Server Badge](https://img.shields.io/github/last-commit/jodyphelan/TBProfiler.svg)](https://github.com/jodyphelan/TBProfiler) [![Codacy Badge](https://app.codacy.com/project/badge/Grade/1f1e0523018c4871a3b56583f59beae9)](https://www.codacy.com/gh/jodyphelan/TBProfiler/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=jodyphelan/TBProfiler&amp;utm_campaign=Badge_Grade) [![Codacy Badge](https://app.codacy.com/project/badge/Coverage/1f1e0523018c4871a3b56583f59beae9)](https://www.codacy.com/gh/jodyphelan/TBProfiler/dashboard?utm_source=github.com&utm_medium=referral&utm_content=jodyphelan/TBProfiler&utm_campaign=Badge_Coverage)



Ce répositoire contient une version complète réécrit de [TB-Profiler web](http://tbdr.lshtm.ac.uk), decrite [ici](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0650-x). Il permet le profilage en utilisant une line de commande d’interface et contient certaines fonctionnalités additionnelles comme la capacité de traiter les données minION.

Le pipeline aligne les séquences reads au génome de référence H37Rv en utilisant bowtie2, BWA ou minimap2 puis identifie les variantes en utilisant bcftools. Ces variantes sont donc comparées a d’autres contenu dans la base de donnes de résistance au médicament. Nous prédisons aussi le nombre de séquences reads supportant les variantes associes a la résistance au médicament comme un aperçu de l’heteroresistance \(pas applicable pour les données de minION\)

## Documentation

Cette page a toutes les informations donc tu as besoin pour démarrer, cependant une documentation additionnelle et plus organisée est disponible [ici](https://jodyphelan.gitbook.io/tb-profiler/). Nous avons aussi une version traduite en d’autre langue:  [:brazil:](https://jodyphelan.gitbook.io/tb-profiler/translations/portugues)[:netherlands:](https://jodyphelan.gitbook.io/tb-profiler/translations/nederlands). Bien vouloir nous contacter si vous aimeriez ameliorer la traduction ou ajouter une nouvelle langue !

## Actualisation
TB-profiler est en développement rapide et constance. Si tu planifie d’utiliser ce programme pour tes travaux, bien vouloir vous assurez que vous utilisez la nouvelle version ! Similairement la base de données n’est pas statique et est continuellement améliore donc assurez-vous d’utiliser la nouvelle version. Dans vos travaux, bien vouloir énoncer la version des deux outils (TB-profiler et la base de données) car ils sont développés de manière indépendante.

## Installation

### Conda

Conda peut fonctionner comme un gestionnaire de packages et est disponible [ici](https://docs.conda.io/en/latest/miniconda.html). Si tu as conda rassure toi que bioconda et conda-forge channels est ajouter :

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Ensuite, vous pouvez installer tb-profiler et toutes ses dépendances a partir du le channel bioconda:

#### Linux

```bash
conda install -c bioconda tb-profiler
```

#### macOS

```bash
conda install -c bioconda tb-profiler 
```

Si vous avez un nouveau mac (M1/M2), vous aurez besoin d’ouvrir votre terminal en utilisant Rosetta. Le majeur parti d’outils bioinformatiques sont construit pour fonctionner sur intel chips (x86) et ne fonctionnera pas sur le M1/M2 chips (ARM). Rosetta imitera un environnement x86 pour vous permettre d’utiliser encore vos outils favorite de bioinformatiques. Pour configurer votre nouvel ordinateur, vous devez :

1.	Ouverez l’application Finder puis allez sur ‘’Application’’
2.	Faites un clique droit sur l’application Terminal
3.	Selectionnez Open using Rosetta
4.	Configurez conda et installer tb-profiler comme décrits ci-dessous

### Instalation manuelle

Il est possible d’installer manuellement, cependant nous vous recommandons grandement d’installer à travers conda. Les dépendances suivantes seront nécessaires au moment de l'exécution : _trimmomatic \(&gt;=v0.38\), bwa \(&gt;=v0.7.17\), minimap2 \(&gt;=v2.16\), samtools \(&gt;=v1.12\), bcftools \(&gt;=v1.12\), freebayes \(&gt;=v1.3.5\), tqdm \(&gt;=v4.32.2\), parallel \(&gt;=v20190522\), samclip \(&gt;=v0.4.0\)_ and snpEff \(&gt;=v5.0.0\)_. Le pipeline devrait fonctionner et été tester avec les programmes dont les versions sont dans les parenthèses.

Pour installer la librairie utilisez les code suivants : 
```bash
pip3 install git+https://github.com/jodyphelan/TBProfiler.git
pip3 install git+https://github.com/jodyphelan/pathogen-profiler.git
mkdir `python -c "import sys; print(getattr(sys, 'base_prefix', getattr(sys, 'real_prefix', sys.prefix)));"`
tb-profiler update_tbdb
```

Vous devriez alors être en mesure d'exécuter en utilisant `tb-profiler`

## Usage

Le premier argument indique le type d’analyse à effectuer. Pour le moment, nous ne prenons en charge que l'identification de petites variantes.

### Exemple de démarrage rapide 

Executez le pipeline entier :

```text
tb-profiler profile -1 /path/to/reads_1.fastq.gz -2 /path/to/reads_2.fastq.gz -p prefix
```

Le préfix es important lorsque vous avez besoin d’exécuter plus d’un échantillon. Ceci conservera les fichier BAM,VCF et résultat dans diffèrent dossier respectives. Les résultats sont en format json et text.

#### Exemple d'exécution

```bash
mkdir test_run
cd test_run
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR166/009/ERR1664619/ERR1664619_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR166/009/ERR1664619/ERR1664619_2.fastq.gz
tb-profiler profile -1 ERR1664619_1.fastq.gz -2 ERR1664619_2.fastq.gz -t 4 -p ERR1664619 --txt
cat results/ERR1664619.results.txt
```

#### Analyse du génome entier

Par défaut, `tb-profiler` analyse uniquement les gènes responsables de la résistance au médicament. Il est aussi possible d’effectuer l’identification des variantes le long du génome entier en utilisant l’argument `--call_whole_genome`. Lorsque cette option est activée, il est également possible de comparer votre échantillon à des séries précédentes pour trouver celles qui sont proches de votre nouvel échantillon. Les échantillons proches sont trouvés en utilisant un seuil de distance SNP qui est activé en utilisant l'argument `--snp_dist`. Les echantillons proches sont stockes sous forme de ficher json et le rapport se retrouvera dans le fichier text.bBien vouloir noter que la fonction `--snp_dist` reste expérimentale et peut produire des résultats non expectes.

#### Spoligotypages

Le spoligotypage experimentale peut être effectues par l’ajoute `--spoligotype` à le command. Ceci est faisable avec les fichiers bam et fastq.

#### Exécutez avec les fichiers BAM existants

En utilisant l’option -a, tu peux spécifier d’utiliser un fichier BAM au lieu des fichiers fastq. **Attention!!!**: le fichier BAM doit avoir été crée en utilisant une version du génome comme base de données qui peut être téléchargeable [ici](ftp://ftp.ensemblgenomes.org/pub/release-32/bacteria//fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/dna/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.toplevel.fa.gz). Ce génome a plusieurs numéros d'accès \(ASM19595v2,NC\_000962.3,GCF\_000195955.2, etc...\). Si vous croyez que votre reference est exactement la même séquence \(la longueur devrait être 4411532\) vous pouvez alors créer une base de données avec le même nom de séquence que celui utilisé. Par exemple, si le nom de votre séquence est "NC\_000962.3" yvous pouvez exécuter la commande suivante avec votre fichier fasta de référence:

```bash
tb-profiler update_tbdb --match_ref /path/to/your/reference/fasta
```

### Résumer l'analyse

Les résultats de nombreux essais peuvent être rassemblés dans un seul tableau à l'aide de la commande suivante :

```bash
tb-profiler collate
```

Cela créera automatiquement un certain nombre de fichiers de résultats regroupés à partir de tous les fichiers de résultats individuels dans le répertoire des résultats. Si vous souhaitez générer ce fichier pour un sous-ensemble de essaies, vous pouvez fournir une liste avec les noms des séries en utilisant l'option `--samples`. Le préfixe des fichiers de sortie est _tbprofiler_ par défaut, mais il peut être modifié avec l'option `--prefix`.

Si `tb-profiler` a été exécuté en utilisant l'argument `--snp_dist`, la fonction `collate` génèrera également un graphique de transmission au format json qui peut être visualisé en glissant et en déposant le fichier dans ce [site](https://jodyphelan.github.io/transmission-graph-viewer)

### Écrire vos propres scripts de résumé

La fonction `collate` extrait les mutations de résistance aux médicaments et la lignée, mais vous pouvez vouloir extraire d'autres caractéristiques qui sont présentes dans les fichiers de résultats json individuels. J'ai créé un petit tutoriel sur la façon de le faire [ici](https://jodyphelan.gitbook.io/tb-profiler/writing-a-custom-collate-script).

## Base de donnees de mutations

TB-Profiler est livré avec une base de données par défaut. Le développement de la bibliothèque de mutations est hébergé sur le [tbdb repository](https://github.com/jodyphelan/tbdb). Veuillez consulter ce repositoire si vous souhaitez vous impliquer dans la base de données ou si vous souhaitez modifier et créer la vôtre.

Si vous souhaitez utiliser une base de données modifiée, vous pouvez télécharger le repositoire tbdb, apporter les modifications nécessaires et exécuter le code suivant à partir du dossier du repositoire tbdb :

```text
tb-profiler create_db --prefix <new_library_name>
tb-profiler load_library --prefix <new_library_name>
```

### Bases de données non-H37Rv

Il est possible d'exécuter TB-Profiler sur un autre génome de référence. Bien qu'il n'y ait actuellement aucun outil d'aide pour créer automatiquement les bases de données pour d'autres références, consultez le [tbdb repository](https://github.com/jodyphelan/tbdb) pour en savoir plus sur ce dont vous avez besoin.

## Sous la hotte

Le pipeline recherche les petits variants et les grandes délétions associés à la résistance aux médicaments. Il rapporte également la lignée. Par défaut, il utilise Trimmomatic pour découper les séquences reads, BWA (ou minimap2 pour nanopore) pour s'aligner sur le génome de référence et bcftools (ou GATKv4/freebayes) pour identifier les variants. Voici un exemple de schéma pour les données paried end d'illumina 

![](https://github.com/jodyphelan/jodyphelan.github.io/raw/master/tb-profiler_uml.svg)

## Fichier ITOL

Plusieurs fichiers sont produits par la fonction `tb-profile collate`. Parmi ceux-ci se trouvent plusieurs fichiers de configuration qui peuvent être utilisés avec [iTOL](http://itol.embl.de/) pour annoter les arbres phylogénétiques. Un petit arbre et des fichiers de configuration ont été placés dans le dossier  _example\_data_ . Pour les utiliser, naviguez sur le site Web d'iTOL et téléchargez le fichier _tbprofiler.tree_ en utilisant le bouton de téléchargement de la barre de navigation. Une fois le fichier téléchargé, vous accédez à une visualisation de l'arbre. Pour ajouter l'annotation, cliquez sur le bouton '+' dans le coin inférieur droit et sélectionnez les fichiers de configuration iTOL. Vous devriez maintenant voir une figure similaire à celle ci-dessous. Les annotations suivantes sont incluses :

*   Lignée
*   Classes de résistance aux médicaments \(Sensitive, RR-TB, HR-TB, MDR-TB, Pre-XDR-TB, XDR\)
*   L’identification de résistance aux médicaments pour les médicaments individuels, les cercles remplis représentent la résistance au médicament.

![](https://github.com/jodyphelan/jodyphelan.github.io/raw/master/img/itol_example.png)

## Citation

i vous prévoyez d'utiliser TB-Profiler dans votre travail, veuillez citer [l'article](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0650-x) et inclure la version de `tb-profiler` et la base de données.

## Problèmes

Veuillez les soulever en utilisant la [problème](https://github.com/jodyphelan/TBProfiler/issues) . Veuillez consulter les [directive relatives aux contributions](https://github.com/jodyphelan/TBProfiler/blob/master/docs/CONTRIBUTING.md) et le [code de conduite](https://github.com/jodyphelan/TBProfiler/blob/master/CODE_OF_CONDUCT.md).

## Questions fréquentes
Q. Comment puis-je utiliser tb-profiler pour produire un arbre phylogénétique ?

A. Pour le moment, TB-Profiler ne fournit aucune fonctionnalité pour créer réellement l'arbre. Les fichiers de configuration pour itol sont générés, mais l'arbre doit être créé avec un autre logiciel. Cela est prévu pour les prochaines versions, mais cela pourrait prendre un certain temps avant d'être complètement développé. Je suggère un autre logiciel pour l'étape de reconstruction phylogénétique (tel que snippy + iqtree).

## Plan future

Veuillez consulter la  [feuille de route](https://github.com/jodyphelan/TBProfiler/projects/1) pour voir les plans futurs ainsi que les fonctionnalités actuellement en cours de développement.
