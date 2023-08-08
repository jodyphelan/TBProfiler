from pathogenprofiler import run_cmd
import pathogenprofiler as pp
import json
import os

collate_text = open("example_collate.txt").read()

if not os.path.isdir("scratch"):
    os.mkdir("scratch")
os.chdir("scratch")

if not os.path.isdir("tb-profiler-test-data"):
    run_cmd("git clone https://github.com/jodyphelan/tb-profiler-test-data.git")

por5_dr_variants = [
    ('rpoB', 'p.Ser450Leu'),
    ('fabG1', 'c.-15C>T'),
    ('inhA', 'p.Ile194Thr'),
    ('pncA', 'p.Val125Gly'),
    ('embB', 'p.Met306Val'),
    ('embB', 'p.Met423Thr'),
    ('gid', 'p.Ala80Pro')
]


def test_db():
    run_cmd("tb-profiler update_tbdb")



def illumina_fastq(caller,mapper):
    run_cmd(f"tb-profiler profile -1 tb-profiler-test-data/por5A.reduced_1.fastq.gz -2 tb-profiler-test-data/por5A.reduced_2.fastq.gz --spoligotype --mapper {mapper} --caller {caller} -p por5A_illumina_{mapper}_{caller}_PE -t 4 --txt --csv ")
    results = json.load(open(f"results/por5A_illumina_{mapper}_{caller}_PE.results.json"))
    return results

def illumina_fastq_single(caller,mapper):
    run_cmd(f"tb-profiler profile -1 tb-profiler-test-data/por5A.reduced_1.fastq.gz --spoligotype --mapper {mapper} --caller {caller} -p por5A_illumina_{mapper}_{caller}_SE -t 4 --txt --csv ")
    results = json.load(open(f"results/por5A_illumina_{mapper}_{caller}_SE.results.json"))
    return results

def test_bwa_freebayes():
    results = illumina_fastq("freebayes","bwa")
    assert results["sublin"] == "lineage4.3.4.2"
    assert results["main_lin"] == "lineage4"
    assert [(v["gene"],v["change"]) for v in results["dr_variants"]] == por5_dr_variants

def test_bwa_bcftools():
    results = illumina_fastq("bcftools","bwa")
    assert results["sublin"] == "lineage4.3.4.2"
    assert results["main_lin"] == "lineage4"
    assert [(v["gene"],v["change"]) for v in results["dr_variants"]] == por5_dr_variants

def test_bwa_gatk():
    results = illumina_fastq("gatk","bwa")
    assert results["sublin"] == "lineage4.3.4.2"
    assert results["main_lin"] == "lineage4"
    assert [(v["gene"],v["change"]) for v in results["dr_variants"]] == por5_dr_variants

def test_bwa_pilon():
    results = illumina_fastq("pilon","bwa")
    assert results["sublin"] == "lineage4.3.4.2"
    assert results["main_lin"] == "lineage4"
    assert [(v["gene"],v["change"]) for v in results["dr_variants"]] == por5_dr_variants

def test_bwa_lofreq():
    results = illumina_fastq("lofreq","bwa")
    assert results["sublin"] == "lineage4.3.4.2"
    assert results["main_lin"] == "lineage4"
    assert [(v["gene"],v["change"]) for v in results["dr_variants"] if v["freq"]>0.05] == por5_dr_variants


def test_vcf1():
    run_cmd("tb-profiler profile -v tb-profiler-test-data/por5A1.vcf.gz --prefix por5_vcf --txt --csv")
    results = json.load(open("results/por5_vcf.results.json"))
    assert results["sublin"] == "lineage4.3.4.2"
    assert results["main_lin"] == "lineage4"
    assert len(results["dr_variants"]) == 7
    assert [(v["gene"],v["change"]) for v in results["dr_variants"]] == por5_dr_variants

def test_collate():
    with open("samples.txt","w") as O:
        O.write("\n".join(["por5A_illumina_bwa_freebayes_PE","por5A_illumina_bwa_gatk_PE","por5A_illumina_bwa_bcftools_PE","por5A_illumina_bwa_pilon_PE","por5A_illumina_bwa_lofreq_PE","por5_vcf"]))
    run_cmd("tb-profiler collate --samples samples.txt")
    assert open("tbprofiler.txt").read() == collate_text
    
def test_minimap2_freebayes():
    results = illumina_fastq("freebayes","minimap2")
    assert results["sublin"] == "lineage4.3.4.2"
    assert results["main_lin"] == "lineage4"
    assert [(v["gene"],v["change"]) for v in results["dr_variants"]] == por5_dr_variants

def test_bwa_freebayes_single():
    results = illumina_fastq_single("freebayes","minimap2")
    assert results["sublin"] == "lineage4.3.4.2"
    assert results["main_lin"] == "lineage4"
    assert [(v["gene"],v["change"]) for v in results["dr_variants"]] == por5_dr_variants

def test_nanopore():
    run_cmd("tb-profiler profile -1 tb-profiler-test-data/por5A.nanopore_reduced.fastq.gz --platform nanopore -p por5A_illumina_nanopore -t 4 --af '0.5,0.7' --depth '0,5' --txt --csv ")
    results = json.load(open("results/por5A_illumina_nanopore.results.json"))
    assert results["sublin"] == "lineage4.3.4.2"
    assert results["main_lin"] == "lineage4"
    assert [(v["gene"],v["change"]) for v in results["dr_variants"]] == por5_dr_variants

def test_fasta():
    run_cmd("tb-profiler profile -f tb-profiler-test-data/por5A1.fasta  -p por5A_fasta --txt --csv")
    results = json.load(open("results/por5A_fasta.results.json"))
    assert results["sublin"] == "lineage4.3.4.2"
    assert results["main_lin"] == "lineage4"
    assert [(v["gene"],v["change"]) for v in results["dr_variants"]] == por5_dr_variants


def test_seqs_from_bam():
    assert pp.get_seqs_from_bam("bam/por5A_illumina_nanopore.bam") == ['Chromosome']