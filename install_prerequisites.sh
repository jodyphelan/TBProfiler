mkdir bin

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

wget http://pathogenseq.lshtm.ac.uk/downloads/TBProfilerFiles.tgz
tar -xvf TBProfilerFiles.tgz
echo "linux" > arch.txt
