mkdir bin

#minimap2
cd minimap2 \
 && make \
 && cp minimap2 ../bin \
 && cd ../ \
 && echo minimap2 OK || echo minimap2 FAIL

#htsbox
cd htsbox \
 && make \
 && mv htsbox ../bin \
 && cd ../ \
 && echo htsbox OK || echo htsbox FAIL

#htslib
cd htslib \
 && autoheader \
 && autoconf \
 && ./configure --disable-lzma;make \
 && mv tabix bgzip ../bin/ \
 && cd ../ \
 && echo htslib OK || echo htslib FAIL

#bwa
cd bwa \
 && make \
 && mv bwa ../bin \
 && cd ../ \
 && echo bwa OK || echo bwa FAIL
#bowtie2
wget -O bowtie2-2.3.4-linux-x86_64.zip https://github.com/BenLangmead/bowtie2/releases/download/v2.3.4/bowtie2-2.3.4-linux-x86_64.zip \
 && unzip bowtie2-2.3.4-linux-x86_64.zip \
 && cd bowtie2-2.3.4-linux-x86_64 \
 && mv bowtie2* ../bin \
 && cd ../ \
 && echo bowtie2 OK || echo bowtie2 FAIL
#bcftools
cd bcftools \
 && make \
 && mv bcftools ../bin \
 && cd ../ \
 && echo bcftools OK || echo bcftools FAIL
#samtools
cd samtools \
 && make \
 && mv samtools ../bin \
 && cd ../ \
 && echo samtools OK || echo samtools FAIL

#lofreq
cd lofreq/dist/ \
 && tar -xvf lofreq_star-2.1.3.1_linux-x86-64.tgz \
 && mv lofreq_star-2.1.3.1/bin/lofreq ../../bin/ \
 && cd ../../ \
 && echo lofreq OK || echo lofreq FAIL

wget http://pathogenseq.lshtm.ac.uk/downloads/TBProfilerFiles.tgz \
 && tar -xvf TBProfilerFiles.tgz \
 && ./bin/bowtie2-build ref/MTB-h37rv_asm19595v2-eg18.fa ref/MTB-h37rv_asm19595v2-eg18.fa
 && echo Files OK || echo Files FAIL
python scripts/create_config.py
