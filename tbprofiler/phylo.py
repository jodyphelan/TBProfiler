import logging
import filelock
import os
from tqdm import tqdm
import pysam
from joblib import Parallel, delayed
import subprocess as sp
from uuid import uuid4
import sys

def run_cmd(cmd: str,log: str=None, desc=None) -> sp.CompletedProcess:
    if desc:
        sys.stderr.write(desc+"\n")
    output = open(log,"w") if log else sp.PIPE
    result = sp.run(cmd,shell=True,check=True,stderr=output,stdout=output)
    return result

def cmd_out(cmd: str) -> str:
    filename = str(uuid4())
    cmd = f"{cmd} > {filename}"
    res = run_cmd(cmd)
    if res.returncode != 0:
        raise Exception(f"Error running {cmd}")
    for line in open(filename):
        yield line.strip()
    os.remove(filename)

def usher_add_sample(args):
    args.bedmask = args.conf['bedmask']

    if args.vcf:
        args.wg_vcf = args.vcf
    else:
        args.wg_vcf = args.files_prefix + ".vcf.gz"
    
    args.tmp_masked_vcf = f"{args.files_prefix}.masked.vcf.gz"
    args.input_phylo = f"{args.dir}/results/phylo.pb"
    args.tmp_output_phylo = f"{args.files_prefix}.pb"
    args.output_nwk = f"{args.files_prefix}.nwk"
    lock = filelock.FileLock(args.input_phylo + ".lock")

    cwd = os.getcwd()
    with lock:
        os.chdir(args.temp)

        tmp_aln = f"{args.files_prefix}.aln.fa"
        tmp_vcf = f"{args.files_prefix}.tmp.vcf"
        consensus_file = prepare_sample_consensus(args.prefix,args)
        run_cmd(f"cat {args.conf['ref']} {consensus_file}> {tmp_aln}")
        run_cmd(f"faToVcf -includeNoAltN {tmp_aln} {tmp_vcf}")
        run_cmd(f"bcftools view -T ^{args.conf['bedmask']} {tmp_vcf} -Oz -o {args.tmp_masked_vcf}")
        run_cmd("usher --vcf %(tmp_masked_vcf)s --load-mutation-annotated-tree %(input_phylo)s --save-mutation-annotated-tree %(tmp_output_phylo)s --write-uncondensed-final-tree" % vars(args))
        run_cmd("mv uncondensed-final-tree.nh %(output_nwk)s" % vars(args))
        run_cmd("matUtils extract -i %(tmp_output_phylo)s -t phylo.nwk" % vars(args))
        for f in ["mutation-paths.txt","placement_stats.tsv"]:
            if os.path.exists(f):
                os.remove(f)
        run_cmd("mv %(tmp_output_phylo)s %(input_phylo)s " % vars(args))
        os.chdir(cwd)

def generate_low_dp_mask(bam,ref,outfile,min_dp = 10):
    missing_positions = []
    ok_positions = set()
    for l in cmd_out(f"samtools depth {bam}"):
        row = l.strip().split("\t")
        if int(row[2])>=min_dp:
            ok_positions.add((row[0],int(row[1])))
    refseq = pysam.FastaFile(ref)
    for chrom,length in zip(refseq.references, refseq.lengths):
        for p in range(length):
            if (chrom,p) not in ok_positions:
                missing_positions.append((chrom,p))
    
    # write missing positions to bed file
    with open(outfile,"w") as O:
        for x in missing_positions:
            O.write(f"{x[0]}\t{x[1]}\t{x[1]+1}\n")

def prepare_usher(treefile,vcf_file):
    run_cmd(f"usher --tree {treefile} --vcf {vcf_file} --collapse-tree --save-mutation-annotated-tree phylo.pb")
    
def prepare_sample_consensus(sample,args):
    s = sample
    tmp_vcf = f"{args.files_prefix}.{s}.vcf.gz"
    run_cmd(f"bcftools norm -m - {args.dir}/vcf/{s}.vcf.gz | bcftools view -T ^{args.bedmask} | bcftools filter --SnpGap 50 | annotate_maaf.py | bcftools filter -S . -e 'MAAF<0.7' |bcftools filter -S . -e 'FMT/DP<20' | bcftools view -v snps -Oz -o {tmp_vcf}")
    run_cmd(f"bcftools index {tmp_vcf}")

    mask_bed = f"{args.files_prefix}.{s}.mask.bed"
    generate_low_dp_mask(f"{args.dir}/bam/{s}.bam",args.conf['ref'],mask_bed)
    run_cmd(f"bcftools consensus --sample {s} -m {mask_bed} -M N -f {args.conf['ref']} {tmp_vcf} | sed 's/>/>{s} /' > {args.files_prefix}.{s}.consensus.fa")
    return f"{args.files_prefix}.{s}.consensus.fa"


def calculate_phylogeny(args):
    samples = [l.strip() for l in open(args.samples)]
    args.tmp_masked_vcf = f"{args.files_prefix}.masked.vcf.gz"
    args.bedmask = args.conf['bedmask']
    
    alignment_file = f"{args.files_prefix}.aln"
    consensus_files = [r for r in tqdm(Parallel(n_jobs=args.threads,return_as='generator')(delayed(prepare_sample_consensus)(s,args) for s in samples),desc="Generating consensus sequences",total=len(samples))]

    run_cmd(f"cat {' '.join(consensus_files)} > {alignment_file}")
    alignment_file_plus_ref = f"{args.files_prefix}.aln.plus_ref"
    run_cmd(f"cat {args.conf['ref']} > {alignment_file_plus_ref}")
    run_cmd(f"cat {alignment_file} >> {alignment_file_plus_ref}")
    tmp_vcf = f"{args.files_prefix}.vcf"
    run_cmd(f"faToVcf {alignment_file_plus_ref} {tmp_vcf}")
    run_cmd(f"iqtree -s {alignment_file} -m GTR+G -nt AUTO",desc="Running IQTree")
    prepare_usher(f"{alignment_file}.treefile",tmp_vcf)
    run_cmd(f"mv phylo.pb {args.dir}/results/")
    os.remove("condensed-tree.nh")
    
        