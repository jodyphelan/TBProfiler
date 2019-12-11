from .utils import filecheck, log, run_cmd
from .bam import bam
from .barcode import barcode, db_compare
from .fastq import fastq
from .abi import abi
from .vcf import vcf
from .fasta import fasta
import re



def bam_profiler(conf, bam_file, prefix, platform, caller, threads=4, whole_genome=False, no_flagstat=False, run_delly=True):

    log("Using %s\n\nPlease ensure that this BAM was made using the same reference as in the database.\nIf you are not sure what reference was used it is best to remap the reads." % bam_file)

    # Put user specified arguments to lower case
    platform = platform.lower()
    caller = caller.lower()

    # Set caller to bcftools if platform is nanopre and wrong caller has been used
    if platform == "nanopore":
        run_delly = False
        caller = "bcftools"

    # Create bam object and call variants
    bam_obj = bam(bam_file, prefix, platform=platform)
    vcf_obj = bam_obj.call_variants(conf["ref"], caller=caller, bed_file=conf["bed"], threads=threads)
    csq = vcf_obj.load_csq(ann_file=conf["ann"])

    # Get % and num reads mapping
    if no_flagstat:
        bam_obj.pct_reads_mapped = "NA"
        bam_obj.num_reads_mapped = "NA"
    else:
        bam_obj.flagstat()


    results = {"variants":[],"qc":{"pct_reads_mapped":bam_obj.pct_reads_mapped,"num_reads_mapped":bam_obj.num_reads_mapped}}
    for sample in csq:
        results["variants"]  = csq[sample]

    mutations = bam_obj.get_bed_gt(conf["barcode"],conf["ref"], caller=caller)
    results["barcode"] = barcode(mutations,conf["barcode"])
    results["missing_regions"] = bam_obj.bed_zero_cov_regions(conf["bed"])

    if run_delly:
        delly_bcf = bam_obj.run_delly()
        deletions = delly_bcf.overlap_bed(conf["bed"])
        for deletion in deletions:
            tmp = {
                "genome_pos": deletion["start"], "gene_id": deletion["region"],
                "chr": deletion["chr"], "freq": 1, "type": "large_deletion",
                "change": "%(chr)s_%(start)s_%(end)s" % deletion
                }
            results["variants"].append(tmp)

    results = db_compare(db_file=conf["json_db"], mutations=results)

    return results


def fasta_profiler(conf, prefix, filename):
    fasta_obj = fasta(filename)
    wg_vcf_file = fasta_obj.get_ref_variants(conf["ref"], prefix, gff=conf["gff"])
    wg_vcf_obj = vcf(wg_vcf_file)
    targets_vcf_file = prefix+".targets.vcf.gz"
    run_cmd("bcftools view -c 1 %s -Oz -o %s -T %s" % (wg_vcf_file,targets_vcf_file,conf['bed']))
    targets_vcf_obj = vcf(targets_vcf_file)
    csq = targets_vcf_obj.load_csq(ann_file=conf["ann"])
    results = {"variants":[],"missing_pos":[],"qc":{"pct_reads_mapped":"NA","num_reads_mapped":"NA"}}
    for sample in csq:
        results["variants"]  = csq[sample]
    mutations = wg_vcf_obj.get_bed_gt(conf["barcode"], conf["ref"])
    barcode_mutations = barcode(mutations,conf["barcode"])
    results["barcode"] = barcode_mutations
    results = db_compare(db_file=conf["json_db"], mutations=results)
    return results


def abi_profiler(conf,prefix,filenames):
    files = filenames.split(",")
    for f in conf:
        filecheck(conf[f])
    for f in files:
        filecheck(f)
    abi_obj = abi(files,prefix)
    vcf_obj = abi_obj.get_variants_vcf(conf["ref"],conf["gff"])
    csq = vcf_obj.load_csq(ann_file=conf["ann"])
    results = {"variants":[]}
    for sample in csq:
        results["variants"]  = csq[sample]
    return results
