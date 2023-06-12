from pathogenprofiler import run_cmd, cmd_out, debug
import filelock
import os
from tqdm import tqdm

def usher_add_sample(args):

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

        run_cmd("bcftools norm -m - %(wg_vcf)s | setGT.py | bcftools view -v snps -Oz -o %(tmp_masked_vcf)s" % vars(args))
        run_cmd("bcftools index %(tmp_masked_vcf)s" % vars(args))
        tmp_aln = f"{args.files_prefix}.aln.fa"
        run_cmd(f"cat {args.conf['ref']} > {tmp_aln}")
        mask_bed = f"{args.files_prefix}.low_dp.mask.bed"
        generate_low_dp_mask(args.bam_file,mask_bed,args.min_depth)
        run_cmd(f"bcftools consensus -m {mask_bed} -M N -f {args.conf['ref']} {args.tmp_masked_vcf} | sed 's/>/>{args.prefix} /' >> {tmp_aln}")
        tmp_vcf = f"{args.files_prefix}.tmp.vcf"
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

def generate_low_dp_mask(bam,outfile,min_dp = 10):
    missing_positions = []
    for l in cmd_out(f"bedtools genomecov -d -ibam {bam}"):
        row = l.strip().split("\t")
        if int(row[2])<min_dp:
            missing_positions.append((row[0],int(row[1])))
    
    # write missing positions to bed file
    with open(outfile,"w") as O:
        for x in missing_positions:
            O.write(f"{x[0]}\t{x[1]}\t{x[1]+1}\n")

def calculate_phylogeny(args):
    samples = [l.strip() for l in open(args.samples)]
    args.tmp_masked_vcf = f"{args.files_prefix}.masked.vcf.gz"
    args.bedmask = args.conf['bedmask']
    
    alignment_file = f"{args.files_prefix}.aln"
    for s in samples:
        tmp_vcf = f"{args.files_prefix}.{s}.vcf.gz"
        run_cmd(f"bcftools norm -m - {args.dir}/vcf/{s}.vcf.gz | bcftools view -T ^{args.bedmask} | setGT.py | bcftools view -v snps -Oz -o {tmp_vcf}")
        run_cmd(f"bcftools index {tmp_vcf}")

        mask_bed = f"{args.files_prefix}.{s}.mask.bed"
        generate_low_dp_mask(f"{args.dir}/bam/{s}.bam",mask_bed)
        run_cmd(f"bcftools consensus -m {mask_bed} -M N -f {args.conf['ref']} {tmp_vcf} | sed 's/>/>{s} /' >> {alignment_file}")
        run_cmd(f"rm {tmp_vcf} {tmp_vcf}.csi")
    alignment_file_plus_ref = f"{args.files_prefix}.aln.plus_ref"
    run_cmd(f"cat {args.conf['ref']} > {alignment_file_plus_ref}")
    run_cmd(f"cat {alignment_file} >> {alignment_file_plus_ref}")
    tmp_vcf = f"{args.files_prefix}.vcf"
    run_cmd(f"faToVcf {alignment_file_plus_ref} {tmp_vcf}")
    run_cmd(f"iqtree -s {alignment_file} -m GTR+G -nt AUTO")
    run_cmd(f"usher --tree {alignment_file}.treefile --vcf {tmp_vcf} --collapse-tree --save-mutation-annotated-tree phylo.pb")

        
    
        