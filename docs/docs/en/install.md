# Install 

TB-Profiler can be installed using a variety of different methods. The easiest is to use conda but you can also install manually.

## Conda 

Conda can function as a package manages are is available [here](https://docs.conda.io/en/latest/miniconda.html). If you have conda make sure the bioconda and conda-forge channels are added:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Then you can install with:

```
conda install -c bioconda tb-profiler
```

If this is taking a long time or complains about inconsistencies there are a few things we can do. The first is switching to `mamba install`. `mamba` basically is a drop-in replacement for conda, but is a lot faster. If you don't have it yet, you can install with:

```
conda install -c conda-forge mamba
```

Second, `tb-profiler` comes bundled with a lot of software. This often can cause compatibility issues with other software installed in your conda base environemnt. It is therefor recommended to intall `tb-profiler` in its own enviroment with:

```
mamba create -n tb-profiler tb-profiler
```

## Manually 

It is possible to install manually. The following pre-requisites will be needed at runtime: trimmomatic (>=v0.38), bwa (>=v0.7.17), minimap2 (>=v2.16), bowtie2 (>=v2.3.5), samtools (>=v1.9), bcftools (>=v1.9), GATK (>=v4.1.4.0), tqdm (>=v4.32.2) and parallel (>=v20190522). The pipeline should work and has been tested on the program versions indicated in parentheses.

To install the library, you can then use:

```
pip3 install git+https://github.com/jodyphelan/TBProfiler.git
pip3 install git+https://github.com/jodyphelan/pathogen-profiler.git
mkdir `python -c "import sys; print(getattr(sys, 'base_prefix', getattr(sys, 'real_prefix', sys.prefix)));"`
tb-profiler update_tbdb
```

## Docker
The pipline can also be run through docker. This is handy when the above methods fail for whatever reason. Pull the latest image using:

```
docker pull quay.io/biocontainers/tb-profiler:4.3.0--pypyh5e36f6f_0
```

If the above method don't work please feel free to [raise an issue](https://github.com/jodyphelan/TBProfiler/issues). We'll try help you get up and running!