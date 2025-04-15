from pathogenprofiler.utils import run_cmd, cmd_out, TempFilePrefix
import pysam
import argparse
import os
from uuid import uuid4
import numpy as np

def robust_bounds(data, k=3):
    median = np.median(data)
    mad = np.median(np.abs(data - median))
    mad_scaled = 1.4826 * mad
    lower = median - k * mad_scaled
    upper = median + k * mad_scaled
    return lower, upper

def generate_low_dp_mask(bam: str,ref: str,outfile: str,min_dp: int = 10) -> None:
    refseq = pysam.FastaFile(ref)
    seqname = refseq.references[0]
    dp = np.zeros(refseq.get_reference_length(seqname))
    for l in cmd_out(f"samtools depth {bam}"):
        row = l.strip().split("\t")
        dp[int(row[1])-1] = int(row[2])

    lower, upper = robust_bounds(dp)
    if min_dp>lower:
        lower = min_dp


    masked_positions = []
    for i, d in enumerate(dp):
        if d < lower or d > upper:
            masked_positions.append(i)

    with open(outfile,"w") as O:
        for i in masked_positions:
            O.write(f"{seqname}\t{i}\t{i+1}\n")



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

def prepare_sample_consensus(
        sample_name: str,
        ref: str,
        input_vcf: str, 
        output_file: str,
        excluded_regions: str,
        low_dp_regions: str = None,

    ) -> str:
    with TempFilePrefix() as tmp:
        tmp_vcf = f"{tmp}.{sample_name}.vcf.gz"
        run_cmd(f"""
            bcftools norm -m - {input_vcf} \
                | bcftools view -T ^{excluded_regions} \
                | annotate_maaf.py \
                | bcftools view -e 'type="indel" && MAAF<0.5' \
                | bcftools filter -S . -e 'GT="alt" && MAAF<0.7' \
                | bcftools filter --SnpGap 50 \
                | rename_vcf_sample.py --sample-name {sample_name} \
                | bcftools view -v snps -Oz -o {tmp_vcf}
        """)
        run_cmd(f"bcftools index {tmp_vcf}")
        if low_dp_regions:
            mask_cmd = f"-m {low_dp_regions} -M N"
        else:
            mask_cmd = ""
        run_cmd(f"bcftools consensus --sample {sample_name} {mask_cmd} -f {ref} {tmp_vcf} | sed 's/>/>{sample_name} /' > {output_file}")

def cli_prepare_sample_consensus(sample: str,input_vcf: str,args: argparse.Namespace) -> str:
    
    mask_bed = f"{args.files_prefix}.{sample}.mask.bed"
    if hasattr(args,'supplementary_bam') and args.supplementary_bam:
        args.bam = args.supplementary_bam
    if args.bam:
        generate_low_dp_mask(f"{args.bam}",args.conf['ref'],mask_bed)
    elif args.low_dp_mask:
        mask_bed = args.low_dp_mask
    elif args.vcf:
        generate_low_dp_mask_vcf(args.vcf,mask_bed)
    else:
        mask_bed = None

    output_file = f'{args.files_prefix}.consensus.fa'

    prepare_sample_consensus(
        sample_name=sample,
        ref=args.conf['ref'],
        input_vcf=input_vcf,
        output_file=output_file,
        excluded_regions=args.conf['bedmask'],
        low_dp_regions=mask_bed,
    )
    return output_file

def get_consensus_vcf(sample: str,input_vcf: str,args: argparse.Namespace) -> str:
    consensus_file = cli_prepare_sample_consensus(sample,input_vcf,args)
    tmp_aln = str(uuid4())
    run_cmd(f"cat {args.conf['ref']} {consensus_file}> {tmp_aln}")
    outfile = f"{args.files_prefix}.masked.vcf"
    run_cmd(f"fa2vcf.py {tmp_aln} {outfile}")
    run_cmd(f'rm {tmp_aln} {tmp_aln}.fai')
    
    return outfile

