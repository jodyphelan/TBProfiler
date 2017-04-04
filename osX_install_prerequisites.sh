mkdir bin
cd htsbox
make
mv htsbox ../bin/
cd ../
cd htslib
make
mv tabix bgzip ../bin/
cd ../
cd bwa
make
mv bwa ../bin
cd ../
wget https://github.com/lomereiter/sambamba/releases/download/v0.6.6/sambamba_v0.6.6_osx.tar.bz2 --no-check-certificate
tar -xvf sambamba_v0.6.6_osx.tar.bz2 
mv sambamba_v0.6.6 bin/sambamba
wget http://pathogenseq.lshtm.ac.uk/downloads/TBProfilerFiles.tgz
tar -xvf TBProfilerFiles.tgz
./bin/snap-aligner index ref/MTB-h37rv_asm19595v2-eg18.fa ref/
echo "mac">arch.txt
