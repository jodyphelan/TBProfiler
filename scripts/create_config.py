import json
import os.path

script_dir = "/".join(os.path.dirname(os.path.realpath(__file__)).split("/")[:-1])
outfile = script_dir+"/conf.json"
conf = {}
conf["annfile"] = "%s/ref/MTB-h37rv_asm19595v2-eg18.tab.ann.gz" % script_dir
conf["bcftools"] = "%s/bin/bcftools" % script_dir
conf["bwa"] = "%s/bin/bwa" % script_dir
conf["bowtie2"] = "%s/bin/bowtie2" % script_dir
conf["samtools"] = "%s/bin/samtools" % script_dir
conf["lofreq"] = "%s/bin/lofreq" % script_dir
conf["htsbox"] = "%s/bin/htsbox" % script_dir
conf["minimap2"] = "%s/bin/minimap2" % script_dir
conf["tabix"] = "%s/bin/tabix" % script_dir
conf["dr_bed_file"] = "%s/db/drdb.bed" % script_dir
conf["dr_json"] = "%s/db/drdb.json" % script_dir
conf["reffile"] = "%s/ref/MTB-h37rv_asm19595v2-eg18.fa" % script_dir
conf["gfffile"] = "%s/ref/MTB-h37rv_asm19595v2-eg18.gff" % script_dir
conf["lineage_bed"] = "%s/db/lineages.bed" % script_dir

for x in conf:
    print(conf[x])
    if not os.path.isfile(conf[x]):
        print("Can't find %s" % x)
        quit()
json.dump(conf,open(outfile,"w"))
