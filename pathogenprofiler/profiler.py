from .utils import filecheck, log, run_cmd
from .bam import bam
from .barcode import barcode, db_compare
from .fastq import fastq
from .vcf import vcf,delly_bcf
from .fasta import fasta
import re



def bam_profiler(conf, bam_file, prefix, platform, caller, threads=1, no_flagstat=False, run_delly=True, calling_params=None, delly_bcf_file=None, run_coverage=True, coverage_fraction_threshold=0, missing_cov_threshold=10, samclip=False, min_depth=10, variant_annotations = False):

    log("Using %s\n\nPlease ensure that this BAM was made using the same reference as in the database.\nIf you are not sure what reference was used it is best to remap the reads." % bam_file)

    ### Put user specified arguments to lower case ###
    platform = platform.lower()
    caller = caller.lower()

    ### Set caller to bcftools if platform is nanopre and wrong caller has been used ###
    if platform == "nanopore":
        run_delly = False
        caller = "bcftools"

    ### Create bam object and call variants ###
    bam_obj = bam(bam_file, prefix, platform=platform, threads=threads)
    vcf_obj = bam_obj.call_variants(conf["ref"], caller=caller, bed_file=conf["bed"], threads=threads, calling_params=calling_params, samclip = samclip, min_dp=missing_cov_threshold)
    if variant_annotations:
        ann_vcf_obj = vcf_obj.add_annotations(conf["ref"],bam_obj.bam_file)
        csq_vcf_obj = ann_vcf_obj.csq(conf["ref"],conf["gff"])
    else:
        csq_vcf_obj = vcf_obj.csq(conf["ref"],conf["gff"])
    csq = csq_vcf_obj.load_csq(ann_file=conf["ann"],extract_ann=True)

    ### Get % and num reads mapping ###
    if no_flagstat:
        bam_obj.pct_reads_mapped = "NA"
        bam_obj.num_reads_mapped = "NA"
    else:
        bam_obj.flagstat()

    ### Put results into a dictionary ###
    results = {
        "variants":[],
        "qc":{
            "pct_reads_mapped":bam_obj.pct_reads_mapped,
            "num_reads_mapped":bam_obj.num_reads_mapped
        }
    }

    if run_coverage:
        results["qc"]["gene_coverage"] = bam_obj.get_region_coverage(conf["bed"], fraction_threshold= coverage_fraction_threshold)
        results["qc"]["missing_positions"] = bam_obj.get_missing_genomic_positions(cutoff=missing_cov_threshold)
    for sample in csq:
        results["variants"]  = csq[sample]

    mutations = bam_obj.get_bed_gt(conf["barcode"],conf["ref"], caller=caller)
    results["barcode"] = barcode(mutations,conf["barcode"])

    ### Run delly if specified ###
    if run_delly:
        if delly_bcf_file:
            delly_bcf_obj = delly_bcf(delly_bcf_file)
        else:
            delly_bcf_obj = bam_obj.run_delly()
        if delly_bcf_obj is not None:
            results["delly"] = "success"
            deletions = delly_bcf_obj.overlap_bed(conf["bed"])
            for deletion in deletions:
                tmp = {
                    "genome_pos": deletion["start"], "gene_id": deletion["region"],
                    "chr": deletion["chr"], "freq": 1, "type": "large_deletion",
                    "change": "%(chr)s_%(start)s_%(end)s" % deletion
                    }
                results["variants"].append(tmp)

        else:
            results["delly"] = "fail"

    ### Compare variants to database ###
    results = db_compare(db_file=conf["json_db"], mutations=results)

    return results


def fasta_profiler(conf, prefix, filename):
    fasta_obj = fasta(filename)
    wg_vcf_file = fasta_obj.get_ref_variants(conf["ref"], prefix)
    wg_vcf_obj = vcf(wg_vcf_file)
    vcf_file = prefix+".targets.vcf.gz"
    run_cmd("bcftools view -c 1 %s -Oz -o %s -T %s" % (wg_vcf_file,vcf_file,conf['bed']))
    vcf_csq_obj = vcf(vcf_file).csq(conf["ref"],conf["gff"])
    csq = vcf_csq_obj.load_csq(ann_file=conf["ann"])
    results = {"variants":[],"missing_pos":[],"qc":{"pct_reads_mapped":"NA","num_reads_mapped":"NA"}}
    for sample in csq:
        results["variants"]  = csq[sample]
    mutations = wg_vcf_obj.get_bed_gt(conf["barcode"], conf["ref"])
    if "C" in mutations["Chromosome"][325505] and  mutations["Chromosome"][325505]["C"]==50:  mutations["Chromosome"][325505] = {"T":25}
    if "G" in mutations["Chromosome"][599868] and  mutations["Chromosome"][599868]["G"]==50:  mutations["Chromosome"][599868] = {"A":25}
    if "C" in mutations["Chromosome"][931123] and  mutations["Chromosome"][931123]["C"]==50:  mutations["Chromosome"][931123] = {"T":25}
    if "T" in mutations["Chromosome"][1759252] and  mutations["Chromosome"][1759252]["T"]==50:  mutations["Chromosome"][1759252] = {"G":25}
    results["barcode"] = barcode(mutations,conf["barcode"])
    results = db_compare(db_file=conf["json_db"], mutations=results)
    run_cmd("rm %s" % wg_vcf_file)
    return results

def vcf_profiler(conf, prefix, sample_name, vcf_file):
    vcf_targets_file = "%s.targets.vcf.gz" % prefix
    run_cmd("bcftools view -T %s %s -Oz -o %s" % (conf["bed"],vcf_file,vcf_targets_file))
    vcf_obj = vcf(vcf_targets_file)
    vcf_csq_obj = vcf_obj.csq(conf["ref"],conf["gff"])
    csq = vcf_csq_obj.load_csq(ann_file=conf["ann"])
    results = {"variants":[],"missing_pos":[],"qc":{"pct_reads_mapped":"NA","num_reads_mapped":"NA"}}
    for sample in csq:
        results["variants"]  = csq[sample]
    mutations = vcf(vcf_file).get_bed_gt(conf["barcode"], conf["ref"])
    if "C" in mutations["Chromosome"][325505] and  mutations["Chromosome"][325505]["C"]==50:  mutations["Chromosome"][325505] = {"T":25}
    if "G" in mutations["Chromosome"][599868] and  mutations["Chromosome"][599868]["G"]==50:  mutations["Chromosome"][599868] = {"A":25}
    if "C" in mutations["Chromosome"][931123] and  mutations["Chromosome"][931123]["C"]==50:  mutations["Chromosome"][931123] = {"T":25}
    if "T" in mutations["Chromosome"][1759252] and  mutations["Chromosome"][1759252]["T"]==50:  mutations["Chromosome"][1759252] = {"G":25}

    results["barcode"] = barcode(mutations,conf["barcode"])
    results = db_compare(db_file=conf["json_db"], mutations=results)
    run_cmd("rm %s* %s*" % (vcf_targets_file,vcf_csq_obj.filename))
    return results
