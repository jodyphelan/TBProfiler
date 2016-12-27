#mkdir bin
#cd htsbox
#make
#mv htsbox ../bin/
#cd ../
#cd htslib
#make
#mv tabix bgzip ../bin/
#cd ../
#cd snap
#make
#mv snap-aligner ../bin/
#cd ../
#wget https://github.com/lomereiter/sambamba/releases/download/v0.6.5/sambamba_v0.6.5_linux.tar.bz2
#tar -xvf sambamba_v0.6.5_linux.tar.bz2
#mv sambamba_v0.6.5 bin/sambamba
wget http://pathogenseq.lshtm.ac.uk/downloads/TBProfilerFiles.tgz
tar -xvf TBProfilerFiles.tgz
