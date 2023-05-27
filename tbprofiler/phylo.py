from pathogenprofiler import run_cmd, errlog
import filelock
import os

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