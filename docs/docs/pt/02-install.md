# Install 

O TB-Profiler pode ser instalado usando uma variedade de métodos diferentes. O mais fácil é por meio de conda, mas você também pode instalar manualmente.

## Conda

Conda pode funcionar como um gerenciador de pacotes e está disponível [aqui](https://docs.conda.io/en/latest/miniconda.html). Se você tiver conda, certifique-se de que os canais bioconda e conda-forge foram adicionados:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Então você pode instalar com:

```
conda install -c bioconda tb-profiler
```

## Manualmente

É possível instalar manualmente. Os seguintes pré-requisitos serão necessários em tempo de execução: trimmomatic (> = v0.38), bwa (> = v0.7.17), minimap2 (> = v2.16), bowtie2 (> = v2.3.5), samtools (> = v1.9), bcftools (> = v1.9), GATK (> = v4.1.4.0), tqdm (> = v4.32.2) e parallel (> = v20190522). O pipeline deve funcionar e já foi testado nas versões do programa indicadas entre parênteses.

Para instalar a biblioteca, você pode usar:

```
pip3 install git+https://github.com/jodyphelan/TBProfiler.git
mkdir `python -c "import sys; print(getattr(sys, 'base_prefix', getattr(sys, 'real_prefix', sys.prefix)));"`
tb-profiler update_tbdb
```