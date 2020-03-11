from __future__ import division
from .bam import bam
from .utils import filecheck, bwa_index, add_arguments_to_self, run_cmd


class fastq:
    """
    Class to hold fastq file and methods.
    Methods include trimming and mapping to a reference genome
    """
    def __init__(self,r1,r2=None,r3=None):
        """
        r1 = Forward reads
        r2 = Reverse reads (optional)
        r3 = Unpaired reads (optional)
        """
        add_arguments_to_self(self, locals())
        # Work out if it is paired end sequencing
        self.paired = True if (r1 and r2) else False
        filecheck(r1)
        self.files = [r1]
        if self.paired:
            filecheck(r2)
            self.files.append(r2)
        if r3:
            filecheck(r3)
            self.files.append(r3)

    def trim(self, prefix, threads=1):
        """Perform trimming"""
        add_arguments_to_self(self, locals())
        if self.paired:
            run_cmd("trimmomatic PE -threads %(threads)s -phred33 %(r1)s %(r2)s -baseout %(prefix)s LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36" % vars(self))
            run_cmd("cat %(prefix)s_1U %(prefix)s_2U > %(prefix)s_TU" % vars(self))
            run_cmd("rm %(prefix)s_1U %(prefix)s_2U" % vars(self))
            return fastq("%(prefix)s_1P" % vars(self), "%(prefix)s_2P" % vars(self), "%(prefix)s_TU" % vars(self))
        else:
            run_cmd("trimmomatic SE -threads %(threads)s -phred33 %(r1)s %(prefix)s_TU LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36" % vars(self))
            return fastq("%(prefix)s_TU" % vars(self))

    def map_to_ref(self, ref_file, prefix, sample_name, aligner, platform, threads=1):
        """Mapping to a reference genome"""
        add_arguments_to_self(self, locals())
        self.aligner = aligner.lower()
        accepted_aligners = ["bwa","bowtie2","minimap2"]
        if self.aligner not in accepted_aligners:
            quit("ERROR: %s not in accepted aligners\n" % aligner)

        self.platform = platform.lower()
        accepted_platforms = ["illumina","nanopore"]
        if self.platform not in accepted_platforms:
            quit("ERROR: %s not in accepted platforms\n" % platform)

        self.bwa_prefix = "bwa mem -t %(threads)s -c 100 -R '@RG\\tID:%(sample_name)s\\tSM:%(sample_name)s\\tPL:%(platform)s' -M -T 50" % vars(self)
        self.bowtie2_prefix = "bowtie2 -p %(threads)s --rg-id '%(sample_name)s' --rg 'SM:%(sample_name)s' --rg 'PL:%(platform)s'" % vars(self)
        self.minimap2_prefix = "minimap2 -t %(threads)s -R '@RG\\tID:%(sample_name)s\\tSM:%(sample_name)s\\tPL:%(platform)s' -a" % vars(self)
        self.bam_file = "%s.bam" % self.prefix
        self.bam_single_file = "%s.single.bam" % self.prefix
        self.bam_pair_file = "%s.pair.bam" % self.prefix
        self.bam_unsort_file = "%s.unsort.bam" % self.prefix
        if self.platform == "nanopore":
            run_cmd("%(minimap2_prefix)s -x map-ont  %(ref_file)s %(r1)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_file)s -" % vars(self))
        else:
            if aligner=="bwa" and self.paired:
                run_cmd("%(bwa_prefix)s %(ref_file)s %(r1)s %(r2)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_pair_file)s -" % vars(self))
                run_cmd("%(bwa_prefix)s %(ref_file)s %(r3)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_single_file)s -" % vars(self))
            elif aligner=="bwa" and not self.paired:
                run_cmd("%(bwa_prefix)s %(ref_file)s %(r1)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_file)s -" % vars(self))
            elif aligner=="bowtie2" and self.paired:
                run_cmd("%(bowtie2_prefix)s -x %(ref_file)s -1 %(r1)s -2 %(r2)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_pair_file)s -" % vars(self))
                run_cmd("%(bowtie2_prefix)s -x %(ref_file)s -U %(r3)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_single_file)s" % vars(self))
            elif aligner=="bowtie2" and not self.paired:
                run_cmd("%(bowtie2_prefix)s  -x %(ref_file)s -1 %(r1)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_file)s -" % vars(self))
            elif aligner=="minimap2" and self.paired:
                run_cmd("%(minimap2_prefix)s -ax sr %(ref_file)s %(r1)s %(r2)s | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_pair_file)s -" % vars(self))
                run_cmd("%(minimap2_prefix)s -ax sr %(ref_file)s %(r3)s| samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_single_file)s -" % vars(self))
            elif aligner=="minimap2" and not self.paired:
                run_cmd("%(minimap2_prefix)s -ax sr %(ref_file)s %(r1)s| samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(bam_file)s -" % vars(self))

            if self.paired:
                run_cmd("samtools merge -@ %(threads)s -f %(bam_unsort_file)s %(bam_pair_file)s %(bam_single_file)s" % vars(self))
                run_cmd("samtools sort -@ %(threads)s -o %(bam_file)s %(bam_unsort_file)s" % vars(self))
                run_cmd("rm %(bam_single_file)s %(bam_pair_file)s %(bam_unsort_file)s" % vars(self))
        return bam(self.bam_file,self.prefix,self.platform)
