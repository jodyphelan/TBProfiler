from .utils import infolog, run_cmd, debug
from .bam import bam
from .barcode import barcode, db_compare
from .vcf import vcf,delly_bcf
from .fasta import fasta




def bam_profiler(conf, bam_file, prefix, platform, caller, threads=1, no_flagstat=False, run_delly=True, calling_params=None, delly_vcf_file=None, run_coverage=True, coverage_fraction_threshold=0, min_depth = 10, missing_cov_threshold=10, samclip=False, variant_annotations = False, call_wg=False):
    infolog("Using %s\n\nPlease ensure that this BAM was made using the same reference as in the database.\nIf you are not sure what reference was used it is best to remap the reads." % bam_file)

    ### Put user specified arguments to lower case ###
    platform = platform.lower()
    caller = caller.lower()

    ### Set caller to bcftools if platform is nanopre and wrong caller has been used ###
    if platform == "nanopore":
        run_delly = False
        caller = "freebayes"

    ### Create bam object and call variants ###
    bam_obj = bam(bam_file, prefix, platform=platform, threads=threads)
    if call_wg:
        wg_vcf_obj = bam_obj.call_variants(conf["ref"], caller=caller, threads=threads, calling_params=calling_params, samclip = samclip, min_dp=min_depth)
        vcf_obj = wg_vcf_obj.view_regions(conf["bed"])
    else:
        vcf_obj = bam_obj.call_variants(conf["ref"], caller=caller, bed_file=conf["bed"], threads=threads, calling_params=calling_params, samclip = samclip, min_dp=min_depth)
    if variant_annotations:
        vcf_obj = vcf_obj.add_annotations(conf["ref"],bam_obj.bam_file)
    else:
        ann_vcf_obj = vcf_obj.run_snpeff(conf["snpEff_db"],conf["ref"],conf["gff"],rename_chroms= conf["chromosome_conversion"])
    ann = ann_vcf_obj.load_ann(bed_file=conf["bed"],upstream=True,synonymous=True,noncoding=True)


    ### Get % and num reads mapping ###
    if no_flagstat:
        bam_obj.pct_reads_mapped = "NA"
        bam_obj.num_reads_mapped = "NA"
        bam_obj.median_coverage = "NA"
    else:
        bam_obj.flagstat()
        bam_obj.get_median_coverage()

    ### Put results into a dictionary ###
    results = {
        "variants":[],
        "qc":{
            "pct_reads_mapped":bam_obj.pct_reads_mapped,
            "num_reads_mapped":bam_obj.num_reads_mapped,
            "median_coverage":bam_obj.median_coverage
        }
    }

    if run_coverage:
        results["qc"]["gene_coverage"] = bam_obj.get_region_coverage(conf["bed"], fraction_threshold= coverage_fraction_threshold)
        results["qc"]["missing_positions"] = bam_obj.get_missing_genomic_positions(cutoff=missing_cov_threshold)
    results["variants"]  = ann

    if "barcode" in conf:
        mutations = bam_obj.get_bed_gt(conf["barcode"],conf["ref"], caller=caller,platform=platform)
        results["barcode"] = barcode(mutations,conf["barcode"])
    

    ### Run delly if specified ###
    if run_delly:
        if delly_vcf_file:
            delly_vcf_obj = delly_bcf(delly_vcf_file)
        else:
            delly_vcf_obj = bam_obj.run_delly()
        if delly_vcf_obj is not None:
            results["delly"] = "success"
            delly_vcf_obj = delly_vcf_obj.get_robust_calls(prefix,conf["bed"])
            ann_vcf_obj = delly_vcf_obj.run_snpeff(conf["snpEff_db"],conf["ref"],conf["gff"],rename_chroms= conf["chromosome_conversion"],split_indels=False)
            results["variants"].extend(ann_vcf_obj.load_ann(bed_file=conf["bed"],ablation=True))
        else:
            results["delly"] = "fail"

    ### Compare variants to database ###
    results = db_compare(db_file=conf["json_db"], mutations=results)
    return results


def fasta_profiler(conf, prefix, filename):
    fasta_obj = fasta(filename)
    wg_vcf_file = fasta_obj.get_ref_variants(conf["ref"], prefix)
    vcf_obj = vcf(wg_vcf_file)
    vcf_file = prefix+".targets.vcf.gz"
    run_cmd("bcftools view -c 1 %s -Oz -o %s -T %s" % (wg_vcf_file,vcf_file,conf['bed']))
    vcf_obj = vcf_obj.run_snpeff(conf["snpEff_db"],conf["ref"],conf["gff"],rename_chroms= conf["chromosome_conversion"])
    ann = vcf_obj.load_ann(bed_file=conf["bed"],upstream=True,synonymous=True,noncoding=True)

    results = {"variants":[],"missing_pos":[],"qc":{"pct_reads_mapped":"NA","num_reads_mapped":"NA"}}
    results["variants"]  = ann
    mutations = vcf_obj.get_bed_gt(conf["barcode"], conf["ref"])
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
    vcf_obj = vcf_obj.run_snpeff(conf["snpEff_db"],conf["ref"],conf["gff"],rename_chroms= conf["chromosome_conversion"])
    ann = vcf_obj.load_ann(bed_file=conf["bed"],upstream=True,synonymous=True,noncoding=True)
    results = {"variants":[],"missing_pos":[],"qc":{"pct_reads_mapped":"NA","num_reads_mapped":"NA"}}
    results["variants"]  = ann
    mutations = vcf(vcf_file).get_bed_gt(conf["barcode"], conf["ref"])
    if "C" in mutations["Chromosome"][325505] and  mutations["Chromosome"][325505]["C"]==50:  mutations["Chromosome"][325505] = {"T":25}
    if "G" in mutations["Chromosome"][599868] and  mutations["Chromosome"][599868]["G"]==50:  mutations["Chromosome"][599868] = {"A":25}
    if "C" in mutations["Chromosome"][931123] and  mutations["Chromosome"][931123]["C"]==50:  mutations["Chromosome"][931123] = {"T":25}
    if "T" in mutations["Chromosome"][1759252] and  mutations["Chromosome"][1759252]["T"]==50:  mutations["Chromosome"][1759252] = {"G":25}

    results["barcode"] = barcode(mutations,conf["barcode"])
    results = db_compare(db_file=conf["json_db"], mutations=results)
    run_cmd("rm %s* %s*" % (vcf_targets_file,vcf_obj.filename))
    return results
