mkdir bin
cd htsbox
make
mv htsbox ../bin/
cd ../
cd htslib
make
mv tabix bgzip ../bin/
cd ../
curl http://snap.cs.berkeley.edu/downloads/snap-1.0beta.1-mac.zip > snap.zip
unzip snap.zip 
mv snap-1.0beta.1-mac/snap bin/snap-aligner
wget https://github.com/lomereiter/sambamba/releases/download/v0.6.6/sambamba_v0.6.6_osx.tar.bz2 --no-check-certificate
tar -xvf sambamba_v0.6.6_osx.tar.bz2 
mv sambamba_v0.6.6 bin/sambamba
wget http://pathogenseq.lshtm.ac.uk/downloads/TBProfilerFiles.tgz
tar -xvf TBProfilerFiles.tgz
