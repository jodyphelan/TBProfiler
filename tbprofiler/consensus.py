from pathogenprofiler.utils import run_cmd,cmd_out
import pysam
import argparse
import os
from uuid import uuid4


def generate_low_dp_mask(bam: str,ref: str,outfile: str,min_dp: int = 10) -> None:
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



def generate_low_dp_mask_vcf(vcf: str,outfile: str,min_dp: int = 10) -> None:
    missing_positions = []
    vcf_obj = pysam.VariantFile(vcf)
    for rec in vcf_obj:
        # use AD field if available
        if 'AD' in rec.samples[0]:
            dp = sum(rec.samples[0]['AD'])
        else:
            dp = rec.samples[0]['DP']
        if dp<min_dp:
            missing_positions.append((rec.chrom,rec.pos))

    # write missing positions to bed file
    with open(outfile,"w") as O:
        for x in missing_positions:
            O.write(f"{x[0]}\t{x[1]}\t{x[1]+1}\n")

def prepare_sample_consensus(sample: str,input_vcf: str,args: argparse.Namespace) -> str:
    s = sample
    tmp_vcf = f"{args.files_prefix}.{s}.vcf.gz"
    run_cmd(f"bcftools norm -m - {input_vcf} | bcftools view -T ^{args.conf['bedmask']} | bcftools view -v snps | annotate_maaf.py | bcftools filter -S . -e 'MAAF<0.7' |bcftools filter -S . -e 'FMT/DP<{args.conf['variant_filters']['depth_soft']}' | bcftools filter --SnpGap 50 |rename_vcf_sample.py --sample-name {s} | bcftools view -v snps -Oz -o {tmp_vcf}")
    run_cmd(f"bcftools index {tmp_vcf}")

    mask_bed = f"{args.files_prefix}.{s}.mask.bed"
    if hasattr(args,'supplementary_bam') and args.supplementary_bam:
        args.bam = args.supplementary_bam
    if args.bam:
        generate_low_dp_mask(f"{args.bam}",args.conf['ref'],mask_bed)
        mask_cmd = f"-m {mask_bed} -M N"
    elif args.low_dp_mask:
        mask_bed = args.low_dp_mask
        mask_cmd = f"-m {mask_bed} -M N"
    elif args.vcf:
        generate_low_dp_mask_vcf(args.vcf,mask_bed)
        mask_cmd = f"-m {mask_bed} -M N"
    else:
        mask_cmd = ""
    run_cmd(f"bcftools consensus --sample {s} {mask_cmd} -f {args.conf['ref']} {tmp_vcf} | sed 's/>/>{s} /' > {args.files_prefix}.{s}.consensus.fa")
    return f"{args.files_prefix}.{s}.consensus.fa"

def get_consensus_vcf(sample: str,input_vcf: str,args: argparse.Namespace) -> str:
    consensus_file = prepare_sample_consensus(sample,input_vcf,args)
    tmp_aln = str(uuid4())
    run_cmd(f"cat {args.conf['ref']} {consensus_file}> {tmp_aln}")
    outfile = f"{args.files_prefix}.masked.vcf"
    run_cmd(f"faToVcf -includeNoAltN {tmp_aln} {outfile}")
    os.remove(tmp_aln)
    return outfile

