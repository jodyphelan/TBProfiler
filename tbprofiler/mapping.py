from files import *
def map(self):
    if self.params["platform"] == "Illumina":
        if self.params["fq1"] and self.params["fq2"]:
            cmd = "%(bwa)s mem -R '@RG\\tID:%(prefix)s\\tSM:%(prefix)s\\tPL:%(platform)s' -t %(threads)s %(reffile)s %(fq1)s %(fq2)s | %(samtools)s view -b -@ %(threads)s - | %(samtools)s sort -@ %(threads)s -o %(bamfile)s -" % self.params
        else:
            cmd = "%(bwa)s mem -R '@RG\\tID:%(prefix)s\\tSM:%(prefix)s\\tPL:%(platform)s' -t %(threads)s %(reffile)s %(fq1)s | %(samtools)s view -b -@ %(threads)s - | %(samtools)s sort -@ %(threads)s -o %(bamfile)s -" % self.params
    elif self.params["platform"] == "minION":
        cmd = "%(minimap2)s -ax map-ont -t %(threads)s %(reffile)s %(fq1)s  | %(samtools)s view -@ %(threads)s -b -  | %(samtools)s sort -@ %(threads)s -o %(bamfile)s - " % self.params
    else: print "Unknown Platform"; quit()
    run_cmd(cmd,verbose=self.params["verbose"])
