# Installation 

TB-Profiler peut être installé en utilisant différentes méthodes. La plus simple est d'utiliser conda mais vous pouvez aussi l'installer manuellement.

## Conda 

Conda peut fonctionner comme un paquet de gestion est disponible [ici] (https://docs.conda.io/en/latest/miniconda.html). Si vous avez conda, assurez-vous que les canaux bioconda et conda-forge sont ajoutés :

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Vous pouvez ensuite l'installer avec :

```
conda install -c bioconda tb-profiler
```

Si cela prend beaucoup de temps ou se plaint d'incohérences, il y a plusieurs choses que nous pouvons faire. La première est de passer à `mamba install`. `mamba` est en fait un remplacement direct de conda, mais il est beaucoup plus rapide. Si vous ne l'avez pas encore, vous pouvez l'installer avec :

```
conda install -c conda-forge mamba
```

Deuxièmement, `tb-profiler` est fourni avec de nombreux logiciels. Cela peut souvent causer des problèmes de compatibilité avec d'autres logiciels installés dans l'environnement de base de Conda. Il est donc recommandé d'installer `tb-profiler` dans son propre environnement avec :

```
mamba create -n tb-profiler tb-profiler
```

## Manuellement 

Il est possible d'effectuer une installation manuelle. Les pré-requis suivants seront nécessaires à l'exécution : trimmomatic (>=v0.38), bwa (>=v0.7.17), minimap2 (>=v2. 16), bowtie2 (>=v2.3.5), samtools (>=v1.9), bcftools (>=v1.9), GATK (>=v4.1.4.0), tqdm (>=v4.32.2) et parallel (>=v20190522). Le pipeline devrait fonctionner et a été testé sur les versions du programme indiquées entre parenthèses.

Pour installer la bibliothèque, vous pouvez ensuite utiliser :

```
pip3 install git+https://github.com/jodyphelan/TBProfiler.git
pip3 install git+https://github.com/jodyphelan/pathogen-profiler.git
mkdir `python -c "import sys; print(getattr(sys, 'base_prefix', getattr(sys, 'real_prefix', sys.prefix)));"`
tb-profiler update_tbdb
```

## Docker
La pipline peut également être exécutée via docker. C'est pratique lorsque les méthodes précédentes échouent pour une raison ou une autre. Téléchargez la dernière image en utilisant :

```
docker pull quay.io/biocontainers/tb-profiler:4.3.0--pypyh5e36f6f_0
```

Si la méthode ci-dessus ne fonctionne pas, n'hésitez pas à [soulever un problème](https://github.com/jodyphelan/TBProfiler/issues). Nous essayerons de vous aider à démarrer !

