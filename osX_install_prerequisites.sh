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
tar -xvf lofreq_star-2.1.3.1_macosx.tgz
mv lofreq_star-2.1.3.1/bin/lofreq ../../bin/
cd ../../

wget http://pathogenseq.lshtm.ac.uk/downloads/TBProfilerFiles.tgz
tar -xvf TBProfilerFiles.tgz
./bin/bwa index ref/MTB-h37rv_asm19595v2-eg18.fa
echo "mac">arch.txt
