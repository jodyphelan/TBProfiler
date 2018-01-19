mkdir bin

#minimap2
cd minimap2
make
cp minimap2 ../bin
cd ../

#htsbox
cd htsbox
make
mv htsbox ../bin
cd ../

#htslib
cd htslib
autoheader
autoconf
./configure --disable-lzma;make
mv tabix bgzip ../bin/
cd ../

#bwa
cd bwa
make
mv bwa ../bin
cd ../

#bowtie2
wget https://github.com/BenLangmead/bowtie2/releases/download/v2.3.4/bowtie2-2.3.4-linux-x86_64.zip
unzip bowtie2-2.3.4-linux-x86_64.zip
cd bowtie2-2.3.4-linux-x86_64
mv bowtie2* ../bin
cd ../

#bcftools
cd bcftools
make
mv bcftools ../bin
cd ../

#samtools
cd samtools
make
mv samtools ../bin
cd ../

#lofreq
cd lofreq/dist/
tar -xvf lofreq_star-2.1.3.1_linux-x86-64.tgz
mv lofreq_star-2.1.3.1/bin/lofreq ../../bin/
cd ../../


wget http://pathogenseq.lshtm.ac.uk/downloads/TBProfilerFiles.tgz
tar -xvf TBProfilerFiles.tgz
./bin/bowtie2-build ref/MTB-h37rv_asm19595v2-eg18.fa ref/MTB-h37rv_asm19595v2-eg18.fa
python scripts/create_config.py
