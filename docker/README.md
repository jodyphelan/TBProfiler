# Docker

To run using docker use the following (assuming you have the data in the current directory).
```
docker pull jodyphelan/tbprofiler
docker run -it -v $PWD:/data jodyphelan/tbprofiler /bin/bash
/opt/TBProfiler/tb-profiler profile -1 sample1_1.fastq.gz -2 sample1_2.fastq.gz
```
