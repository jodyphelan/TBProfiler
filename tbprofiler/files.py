import gzip
import os.path
import subprocess

def filecheck(filename):
    if not os.path.isfile(filename):
        print "Can't find %s" % filename
        quit()
    else:
        return True

def run_cmd(cmd,verbose=False):
    if verbose==2:
        print "\nRunning command:\n%s" % cmd
        stderr = open("/dev/stderr","w")
    elif verbose==1:
        print "\nRunning command:\n%s" % cmd
        stderr = open("/dev/null","w")
    else:
        stderr = open("/dev/null","w")

    res = subprocess.call(cmd,shell=True,stderr = stderr)
    stderr.close()
    if res!=0:
        print "Command Failed! Please Check!"
        quit()

def init_storage(self):
    for x in ["bam","vcf","results"]:
        tdir = self.params["stor_dir"]+"/"+x
        if not os.path.isdir(tdir):
            os.mkdir(tdir)

def index_bam(self):
    if os.path.isfile("%(bamfile)s.bai" % self.params):
        cmd = "%(samtools)s index %(bamfile)s" % self.params
        run_cmd(cmd)

def cleanup(self):
    if self.params["mapping"]:
        cmd = "mv %(bamfile)s* %(stor_dir)s/bam/ && mv %(dr_vcffile)s %(stor_dir)s/vcf/ && mv %(json_results)s %(txt_results)s %(stor_dir)s/results/ && rm %(prefix)s.*" % self.params
    else:
        cmd = "mv %(dr_vcffile)s %(stor_dir)s/vcf/ && mv %(json_results)s %(txt_results)s %(stor_dir)s/results/ && rm %(prefix)s.*" % self.params
    run_cmd(cmd)

def verify_fq(filename):
    FQ = open(filename) if filename[-3:]!=".gz" else gzip.open(filename)
    l1 = FQ.readline()
    if l1[0]!="@":
        print "First character is not \"@\"\nPlease make sure this is fastq format\nExiting..."
        quit()
    else:
        return True
