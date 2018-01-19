from files import *
def map(self):
	if self.params["platform"] == "Illumina":
		if self.mapper == "bwa":
			if self.params["fq1"] and self.params["fq2"]:
				cmd = "%(bwa)s mem -R '@RG\\tID:%(prefix)s\\tSM:%(prefix)s\\tPL:%(platform)s' -t %(threads)s %(reffile)s %(fq1)s %(fq2)s | %(samtools)s view -b -@ %(threads)s - | %(samtools)s sort -@ %(threads)s -o %(bamfile)s -" % self.params
			else:
				cmd = "%(bwa)s mem -R '@RG\\tID:%(prefix)s\\tSM:%(prefix)s\\tPL:%(platform)s' -t %(threads)s %(reffile)s %(fq1)s | %(samtools)s view -b -@ %(threads)s - | %(samtools)s sort -@ %(threads)s -o %(bamfile)s -" % self.params
		elif self.mapper == "bowtie2":
			if self.params["fq1"] and self.params["fq2"]:
				cmd = "%(bowtie2)s -x %(reffile)s -1 %(fq1)s -2 %(fq2)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bamfile)s -" % self.params
			if self.params["fq1"]:
				cmd = "%(bowtie2)s -x %(reffile)s -U %(fq1)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bamfile)s -" % self.params
	elif self.params["platform"] == "minION":
		cmd = "%(minimap2)s -ax map-ont -t %(threads)s %(reffile)s %(fq1)s  | %(samtools)s view -@ %(threads)s -b -  | %(samtools)s sort -@ %(threads)s -o %(bamfile)s - " % self.params
	else:
		print "Unknown Platform"
		quit()
	print cmd
	run_cmd(cmd,verbose=self.params["verbose"])
