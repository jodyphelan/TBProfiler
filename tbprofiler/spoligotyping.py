import os 
import pathogenprofiler as pp

def bam2spoligotype(bamfile,files_prefix,conf):
    chrom = open(conf['bed']).readline().split()[0]
    pp.run_cmd(f"samtools view -b {bamfile} {chrom}:3117003-3127206 | samtools fastq > {files_prefix}.spacers.fq")
    fastq = pp.fastq(files_prefix+".spacers.fq")
    kmers = fastq.get_kmer_counts(files_prefix,klen=25)
    counts = kmers.load_kmer_counts(conf['spoligotype_spacers'])
    binary,octal = counts2spoligotype(counts)
    os.remove(f"{files_prefix}.spacers.fq")
    return {"binary":binary,"octal":octal,"spacers":counts}

def counts2spoligotype(counts,cutoff=None):
    spacers = []
    if cutoff==None:
        cutoff = min([10,max([x["count"] for x in counts])*0.2])
    for k in counts:
        spacers.append("1" if k['count']>cutoff else "0")
    
    octal = []
    for i in range(0,40,3):
        tmp = "".join([spacers[i],spacers[i+1],spacers[i+2]])
        if tmp=="000":octal.append("0")
        elif tmp=="001":octal.append("1")
        elif tmp=="010":octal.append("2")
        elif tmp=="011":octal.append("3")
        elif tmp=="100":octal.append("4")
        elif tmp=="101":octal.append("5")
        elif tmp=="110":octal.append("6")
        elif tmp=="111":octal.append("7")

    octal.append("0" if spacers[42]=="0" else "1")
    sitvit_str = "".join(["n" if x=="1" else "o" for x in spacers])
    binary_str = "".join(spacers)
    octal_str = "".join(octal)
    return binary_str,octal_str
