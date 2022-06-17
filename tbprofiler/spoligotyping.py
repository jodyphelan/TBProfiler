import os
import pathogenprofiler as pp
import csv

def spoligotype(args):
    if "bam_file" in vars(args) and args.bam_file:
        result = bam2spoligotype(args.bam_file,args.files_prefix,args.conf)
    elif args.read1:
        result = fq2spoligotype(args.read1,args.files_prefix,args.conf,args.read2)
    elif args.fasta:
        result = fa2spoligotype(args.fasta,args.files_prefix,args.conf)
    ann = get_spoligotype_annotation(result["octal"],args.conf['spoligotype_annotations'])
    result.update(ann)
    return result

def fa2spoligotype(fasta,files_prefix,conf):
    fasta = pp.fasta(fasta)
    kmers = fasta.get_kmer_counts(files_prefix,klen=25)
    counts = kmers.load_kmer_counts(conf['spoligotype_spacers'])
    binary,octal = counts2spoligotype(counts,cutoff=1)
    return {"binary":binary,"octal":octal,"spacers":counts}


def fq2spoligotype(r1,files_prefix,conf,r2=None):
    fastq = pp.fastq(r1,r2)
    kmers = fastq.get_kmer_counts(files_prefix,klen=25)
    counts = kmers.load_kmer_counts(conf['spoligotype_spacers'])
    binary,octal = counts2spoligotype(counts)
    return {"binary":binary,"octal":octal,"spacers":counts}

def bam2spoligotype(bamfile,files_prefix,conf):
    chrom = open(conf['bed']).readline().split()[0]
    tmp_fq_file = f"{files_prefix}.spacers.fq"
    pp.run_cmd(f"samtools view -b {bamfile} {chrom}:3117003-3127206 | samtools fastq > {tmp_fq_file}")
    results = fq2spoligotype(tmp_fq_file,files_prefix,conf)
    os.remove(tmp_fq_file)
    return results

def get_spoligotype_annotation(octal,annotation_csv):
    result = {
        "family": None,
        "SIT": None,
        "countries": None
    }
    for row in csv.DictReader(open(annotation_csv)):
        if row["Spoligo Octal"][1:]==octal:
            result["family"] = row["Lineage (SITVIT2)"]
            result["SIT"] = row["SIT"]
            result["countries"] = row["Country Distribution (SITVIT2)"]
    return result

def counts2spoligotype(counts,cutoff=None):
    spacers = []
    if cutoff==None:
        cutoff = min([10,max([x["count"] for x in counts])*0.2])
    for k in counts:
        spacers.append("1" if k['count']>=cutoff else "0")
    
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

