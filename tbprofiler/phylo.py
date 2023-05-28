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
    args.bedmask = args.conf['bedmask']
    lock = filelock.FileLock(args.input_phylo + ".lock")

    cwd = os.getcwd()
    with lock:
        run_cmd("bcftools view -T ^%(bedmask)s %(wg_vcf)s -Oz -o %(tmp_masked_vcf)s" % vars(args))
        os.chdir(args.temp)
        run_cmd("/Users/jody/miniconda3/envs/usher/bin/usher --vcf %(tmp_masked_vcf)s --load-mutation-annotated-tree %(input_phylo)s --save-mutation-annotated-tree %(tmp_output_phylo)s --write-uncondensed-final-tree" % vars(args))
        run_cmd("mv uncondensed-final-tree.nh %(output_nwk)s" % vars(args))
        run_cmd("rm mutation-paths.txt placement_stats.tsv")
        run_cmd("mv %(tmp_output_phylo)s %(input_phylo)s " % vars(args))
        os.chdir(cwd)

def get_nucleotide(row,min_cov = 10,min_freq = 0.8):
    '''Get nucleotide from row'''
    pos, sample, gt, tgt, ad = row
    if gt==".": 
        return "N"
    ad = [int(x) for x in ad.split(",")]
    if sum(ad)<=min_cov:
        return "N"
    adf = sorted([float(x/sum(ad)) for x in ad])
    if adf[-1]<min_freq:
        return "N"
    return tgt[0]

def dump_buffer(buffer,args):
    '''Dump buffer to file'''
    variant_positions = []
    for i in range(len(list(buffer.values())[0])):
        if len(set([buffer[s][i] for s in buffer]) - set(["N"]))>1:
            variant_positions.append(i)
    for s in buffer:
        with open(f"{args.files_prefix}.{s}.fa","a") as O:
            O.write("".join([buffer[s][i] for i in variant_positions]))

def get_sample_name_from_bam_header(bam):
    '''Get sample name from bam header'''
    sample_name = None
    for l in cmd_out(f"samtools view -H {bam}"):
        if l.startswith("@RG"):
            sample_name =  l.split("\t")[1].split(":")[1]
    return sample_name

def calculate_phylogeny(args):
    samples = [l.strip() for l in open(args.samples)]
    args.tmp_masked_vcf = f"{args.files_prefix}.masked.vcf.gz"
    args.bedmask = args.conf['bedmask']
    
    alignment_file = f"{args.files_prefix}.aln"
    for s in samples:
        tmp_vcf = f"{args.files_prefix}.{s}.vcf.gz"
        run_cmd(f"bcftools norm -m - {args.dir}/vcf/{s}.vcf.gz | bcftools view -T ^{args.bedmask} | setGT.py | bcftools view -v snps -Oz -o {tmp_vcf}")
        run_cmd(f"bcftools index {tmp_vcf}")
        run_cmd(f"bcftools consensus -f {args.conf['ref']} {tmp_vcf} | sed 's/\>/\>{s} /' >> {alignment_file}")
        run_cmd(f"rm {tmp_vcf} {tmp_vcf}.csi")
    alignment_file_plus_ref = f"{args.files_prefix}.aln.plus_ref"
    run_cmd(f"cat {args.conf['ref']} > {alignment_file_plus_ref}")
    run_cmd(f"cat {alignment_file} >> {alignment_file_plus_ref}")
    tmp_vcf = f"{args.files_prefix}.vcf"
    run_cmd(f"faToVcf {alignment_file_plus_ref} {tmp_vcf}")
    run_cmd(f"iqtree -s {alignment_file} -m GTR+G -nt AUTO")
    run_cmd(f"usher --tree {alignment_file}.treefile --vcf {tmp_vcf} --collapse-tree --save-mutation-annotated-tree phylo.pb")

        
    
        