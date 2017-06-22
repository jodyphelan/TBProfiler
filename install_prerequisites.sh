
#htslib
cd htslib
autoheader
autoconf
./configure --disable-lzma;make
mv tabix bgzip ../bin/
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

#snap
cd snap
make
mv snap-aligner ../bin/
cd ../

#lofreq
cd lofreq/dist/
tar -xvf lofreq_star-2.1.3.1_linux-x86-64.tgz 
mv lofreq_star-2.1.3.1/bin/lofreq ../../bin/
cd ../../


wget http://pathogenseq.lshtm.ac.uk/downloads/TBProfilerFiles.tgz
tar -xvf TBProfilerFiles.tgz
echo "linux" > arch.txt
