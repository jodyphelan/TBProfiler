from pathogenprofiler import run_cmd
import pathogenprofiler as pp
import json
import os
import pytest

collate_text = open("example_collate.txt").read()

if not os.path.isdir("scratch"):
    os.mkdir("scratch")
os.chdir("scratch")

if not os.path.isdir("tb-profiler-test-data"):
    run_cmd("git clone https://github.com/jodyphelan/tb-profiler-test-data.git")


por5_dr_variants = [
    ('rpoB', 'p.Ser450Leu'),
    ('inhA', 'c.-777C>T'),
    ('pncA', 'p.Val125Gly'),
    ('embB', 'p.Met306Val'),
]

db = 'testdb'
branch = 'who'

def test_db():
    run_cmd(f"tb-profiler update_tbdb --branch {branch} --prefix {db}")


def check_assertations(filename):
    results = json.load(open(filename))
    assert results["sub_lineage"] == "lineage4.3.4.2"
    assert results["main_lineage"] == "lineage4"
    assert [(v["gene_name"],v["change"]) for v in results["dr_variants"]] == por5_dr_variants


def test_fastq():
    run_cmd(f"tb-profiler profile --db {db} -1 tb-profiler-test-data/por5A.reduced_1.fastq.gz -2 tb-profiler-test-data/por5A.reduced_2.fastq.gz -p por5A_fastq -t 4 --txt --csv ")
    check_assertations("results/por5A_fastq.results.json")

variant_callers = ['bcftools','freebayes','gatk','lofreq','pilon']
@pytest.mark.parametrize("caller",variant_callers)
def test_variant(caller):
    run_cmd(f"tb-profiler profile --db {db} -a bam/por5A_fastq.bam --caller {caller} -p por5A_{caller} -t 4 --txt --csv --docx ")
    check_assertations(f"results/por5A_{caller}.results.json")

def test_vcf():
    run_cmd(f"tb-profiler profile --db {db} -v tb-profiler-test-data/por5A1.vcf.gz --prefix por5_vcf --txt --csv --docx")
    check_assertations("results/por5_vcf.results.json")
    
def test_nanopore():
    run_cmd(f"tb-profiler profile --db {db} -1 tb-profiler-test-data/por5A.nanopore_reduced.fastq.gz --platform nanopore -p por5A_nanopore -t 4 --caller bcftools --af '0.5,0.7' --depth '0,5' --txt --csv --docx")
    check_assertations("results/por5A_nanopore.results.json")

def test_fasta():
    run_cmd(f"tb-profiler profile --db {db} -f tb-profiler-test-data/por5A1.fasta  -p por5A_fasta --txt --csv --docx")
    results = json.load(open('results/por5A_fasta.results.json'))
    assert results["sub_lineage"] == "lineage4.3.4.2"
    assert results["main_lineage"] == "lineage4"
    assert set([(v["gene_name"],v["change"]) for v in results["dr_variants"]]) - set([('fbiC','c.2565_*117del')]) == set(por5_dr_variants)


def test_collate():
    with open("samples.txt","w") as O:
        for caller in variant_callers:
            O.write(f"por5A_{caller}\n")
        O.write("por5_vcf\n")
    run_cmd(f"tb-profiler collate --db {db} --samples samples.txt")
    assert open("tbprofiler.txt").read() == collate_text

def test_tbp_parser():
    run_cmd("git clone https://github.com/theiagen/tbp-parser.git")
    run_cmd("samtools index bam/por5A_fastq.bam")
    run_cmd('python tbp-parser/tbp_parser/tbp_parser.py results/por5A_freebayes.results.json bam/por5A_fastq.bam    -o "example-tbp-parser"    --min_depth 12     --min_frequency 0.9     --sequencing_method "Illumina NextSeq"    --operator "John Doe"')
    test_text = ','.join([l.strip().split(',')[:-2] for l in open("example-tbp-parser.looker_report.csv")][1])
    target_text = "por5A_freebayes,Illumina NextSeq,U,S-Interim,U,U,R,R,R,U,S,U,S,R,R,U,lineage4,DNA of Mycobacterium tuberculosis species detected"
    assert test_text == target_text