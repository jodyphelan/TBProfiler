mkdir bin

#htslib
cd htslib
autoheader
autoconf
./configure --disable-lzma;make
mv tabix bgzip ../bin/
cd ../

#bwa
wget http://pathogenseq.lshtm.ac.uk/downloads/bwa-0.6.2.tgz
tar -xvf bwa-0.6.2.tgz
cd bwa-0.6.2; make
mv bwa ../bin

#bcftools
cd bcftools
make
mv bcftools ../bin
cd ../

#samtools
git clone https://github.com/samtools/samtools.git
cd samtools
make
mv samtools ../bin
cd ../

wget http://pathogenseq.lshtm.ac.uk/downloads/TBProfilerFiles.tgz
tar -xvf TBProfilerFiles.tgz
./bin/bwa index ref/MTB-h37rv_asm19595v2-eg18.fa
echo "rpi">arch.txt
