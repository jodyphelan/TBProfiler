mkdir bin

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
echo "mac">arch.txt
