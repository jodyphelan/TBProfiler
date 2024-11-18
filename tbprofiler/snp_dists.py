import pickle
from pathogenprofiler import cmd_out
import logging
import os
import json
from .output import write_outputs
from copy import copy
import filelock
import sqlite3
from tqdm import tqdm
from .consensus import get_consensus_vcf
import argparse
from .models import ProfileResult, LinkedSample
from typing import List, Tuple

def extract_variant_set(vcf_file: str) -> Tuple[set,set]:
    ref_diffs = set()
    missing = set()
    for l in cmd_out(f"bcftools view {vcf_file} | bcftools query -f '%POS[\t%GT]\n'"):
        if l[0]=="#": continue
        row = l.strip().split()
        pos = int(row[0])
        gt = row[1]
        if gt==".": 
            missing.add(pos)
            continue
        elif gt=="1":
            ref_diffs.add(int(pos))
        else:
            raise Exception("Unknown GT: %s" % gt)

    return ref_diffs, missing

class DB:
    """
    Class for storing and searching for SNP differences between samples
    
    Arguments
    ---------
    filename : str
        Filename of the sqlite database to use
    
    Attributes
    ----------
    filename : str
        Filename of the sqlite database to use
    conn : sqlite3.Connection
        Connection to the sqlite database
    c : sqlite3.Cursor
        Cursor to the sqlite database
        
    Methods
    -------
    store(result,vcf_file)
        Store the SNP differences from a sample in the database
    search(result,vcf_file,cutoff=20)
        Search for samples with similar SNP differences in the database
    """
    def __init__(self, filename: str) -> None:
        self.filename = filename
        self.conn = sqlite3.connect(filename)
        self.c = self.conn.cursor()
        self.c.execute('''CREATE TABLE IF NOT EXISTS variants (sample text, lineage text, diffs binary, missing binary)''')
    def store(self,result: ProfileResult, vcf_file: str) -> None:
        sample_name = result.id
        diffs,missing = extract_variant_set(vcf_file)
        res = self.c.execute("SELECT sample FROM variants WHERE sample=?",(sample_name,)).fetchone()
        if res:
            self.c.execute("UPDATE variants SET lineage=?, diffs=?, missing=? WHERE sample=?",(result.sub_lineage, pickle.dumps(diffs), pickle.dumps(missing), sample_name))
        else:
            self.c.execute("INSERT INTO variants VALUES (?,?,?,?)",(sample_name, result.sub_lineage, pickle.dumps(diffs), pickle.dumps(missing)))
        self.conn.commit()
        self.diffs = diffs
        self.missing = missing
    def search(self,result: ProfileResult, vcf_file: str, cutoff: int = 20) -> List[LinkedSample]:
        logging.info("Searching for close samples in %s" % self.filename)
        self.c.execute("SELECT sample, diffs, missing FROM variants WHERE lineage=?",(result.sub_lineage,))
        self.diffs,self.missing = extract_variant_set(vcf_file)
        sample_dists = []
        for s,d,m in tqdm(self.c.fetchall(),desc="Searching for close samples"):
            dist = self.diffs.symmetric_difference(pickle.loads(d))
            dist -= self.missing
            dist -= pickle.loads(m) 
            if (ld:=len(dist))<cutoff:
                sample_dists.append(
                    LinkedSample(
                        sample = s,
                        distance = ld,
                        positions = list(dist)
                    )
                )
        logging.info("Found %s close samples" % len(sample_dists))
        return sample_dists

def sample_in_linked_list(sample,result_file):
    result = ProfileResult(**json.load(open(result_file)))
    if sample in [d.sample for d in result.linked_samples]:
        return True
    else:
        return False

def run_snp_dists(args: argparse.Namespace,result: ProfileResult) -> None:
    logging.info("Calculating SNP distances to look for closely related samples")
    if args.vcf:
        wg_vcf = args.vcf
    else:
        wg_vcf = args.files_prefix + ".vcf.gz"
    input_vcf = get_consensus_vcf(args.prefix, wg_vcf,args)
    logging.debug("Input VCF: %s" % input_vcf)
    if args.snp_diff_db:
        dbname = args.snp_diff_db
    else:
        dbname = f'{args.dir}/results/snp_diffs.db'
    lock = f"{dbname}.lock"
    with filelock.SoftFileLock(lock):
        db = DB(dbname)
        linked_samples = db.search(result,input_vcf,args.snp_dist)
        if not args.snp_diff_no_store:
            db.store(result,input_vcf)
        result.linked_samples = [d for d in linked_samples if d.sample!=result.id]



def update_neighbour_snp_dist_output(args: argparse.Namespace,result: ProfileResult) -> None:
    for s in result.linked_samples:
        logging.debug("Close sample found: %s (%s). Updating result files" % (s.sample,s.distance))
        f = os.path.join(args.dir,"results",f"{s.sample}.results.json")
        if not os.path.exists(f):
            continue
        if not sample_in_linked_list(args.prefix,f):
            lock = filelock.SoftFileLock(f + ".lock")
            with lock:
                logging.debug("Acquiring lock for %s" % f)
                data = ProfileResult(**json.load(open(f)))
                data.linked_samples.append(
                    LinkedSample(
                        sample = args.prefix,
                        distance = s.distance,
                        positions = s.positions
                    )
                )
                temp_args = copy(args)
                temp_args.prefix = s.sample
                write_outputs(temp_args,data,template_file=args.text_template)
                logging.debug("Finished with lock for %s" % f)