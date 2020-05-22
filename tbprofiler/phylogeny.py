import json
import pathogenprofiler as pp
import os
def phylogeny(prefix,conf,sample_file=None,base_dir = ".",threads=3):
    if sample_file:
        samples = [x.rstrip() for x in open(sample_file).readlines()]
    else:
        samples = [x.replace(".results.json","") for x in os.listdir("results/") if x[-13:]==".results.json"]

    samples_file = pp.get_random_file()
    OUT = open(samples_file,"w")
    OUT.write("%s\n"%"\n".join(samples))
    OUT.close()
    for s in samples:
        tprefix = s+".genome"
        gbcf_file = "%s.gbcf" % tprefix
        if pp.nofile("%s/vcf/%s.genome.gbcf" % (base_dir,s)):
            bam_file = "%s/bam/%s.bam" % (base_dir,s)
            bam_obj = pp.bam(bam_file,s,conf["ref"])
            bam_obj.gbcf(prefix=tprefix)
            pp.run_cmd("mv %s* %s/vcf" % (gbcf_file,base_dir))
    cmd = "merge_vcfs.py %s %s %s --vcf_dir %s/vcf/ --vcf_ext genome.gbcf" % (samples_file,conf["ref"],prefix,base_dir)
    print(cmd)
    #pp.rm_files([samples_file])
