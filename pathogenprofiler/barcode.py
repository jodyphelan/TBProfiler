from .utils import *
import json
import re

def get_missense_codon(x):
    re_obj = re.search("([0-9]+)",x)
    if re_obj:
        return int(re_obj.group(1))
    else:
        log("Error can't find codon number in %s" % x,True)

def get_indel_nucleotide(x):
    re_obj = re.search("([0-9]+)",x)
    if re_obj:
        return int(re_obj.group(1))
    else:
        log("Error can't find nucleotide number in %s" % x,True)

def barcode(mutations,barcode_bed):
    bed_num_col = len(open(barcode_bed).readline().rstrip().split())
    cols = [1,3,4,5,6]+list(range(7,bed_num_col+1)) if bed_num_col>6 else [1,3,4,5,6]
    bed = []
    lineage_info = {}
    for l in open(barcode_bed):
        row = l.strip().split()
        bed.append(row)
        lineage_info[row[3]] = row
    #{'Chromosome':{'4392120': ('Chromosome', '4392120', 'lineage4.4.1.2', 'G', 'A', 'Euro-American', 'T1', 'None')}}
    barcode_support = defaultdict(list)
    for marker in bed:
        tmp = [0,0]
        chrom,pos = marker[0],int(marker[2])
        if chrom in mutations and pos in mutations[chrom]:
            if marker[4] in mutations[chrom][pos]:
                tmp[1] = mutations[chrom][pos][marker[4]]
            tmp[0] = sum(list(mutations[chrom][pos].values())) - tmp[1]

        if  tmp==[0,0]: continue
        barcode_support[marker[3]].append(tmp)
    barcode_frac = defaultdict(float)
    for l in barcode_support:
        # if stdev of fraction across all barcoding positions > 0.15
        if stdev([x[1]/(x[0]+x[1]) for x in barcode_support[l]])>0.15: continue
        # if number of barcoding positions > 5 and only one shows alternate
        if len(barcode_support[l])>5 and len([x for x in barcode_support[l] if (x[1]/(x[0]+x[1]))>0])<2: continue
        barcode_pos_reads = sum([x[1] for x in barcode_support[l]])
        barcode_neg_reads = sum([x[0] for x in barcode_support[l]])
        lf = barcode_pos_reads/(barcode_pos_reads+barcode_neg_reads)
        if lf<0.05:continue
        barcode_frac[l] = lf
    final_results = []

    for l in barcode_frac:
        tmp = {"annotation":l,"freq":barcode_frac[l],"info":[]}
        if bed_num_col>6:
            tmp["info"] = [lineage_info[l][i] for i in range(5,bed_num_col-1)]
        final_results.append(tmp)
    return final_results

def db_compare(mutations,db_file):
    db = json.load(open(db_file))
    annotated_results = mutations
    for i in range(len(mutations["variants"])):
        #var = {'genome_pos': 6140, 'gene_id': 'Rv0005', 'chr': 'Chromosome', 'freq': 0.975609756097561, 'type': 'missense', 'change': '301V>301L'}
        var = mutations["variants"][i]
        if var["gene_id"] in db:
            db_var_match = None
            if var["change"] in db[var["gene_id"]]:
                db_var_match = db[var["gene_id"]][var["change"]]
            elif "frameshift" in var["type"] and "frameshift" in db[var["gene_id"]]:
                db_var_match = db[var["gene_id"]]["frameshift"]
            elif "missense" in var["type"] and "any_missense_codon_%s" % get_missense_codon(var["change"]) in db[var["gene_id"]]:
                db_var_match = db[var["gene_id"]]["any_missense_codon_%s" % get_missense_codon(var["change"])]
            elif "frame" in var["type"] and "any_indel_nucleotide_%s" % get_indel_nucleotide(var["change"]) in db[var["gene_id"]]:
                db_var_match = db[var["gene_id"]]["any_indel_nucleotide_%s" % get_indel_nucleotide(var["change"])]
            elif "stop_gained" in var["type"] and "premature_stop" in db[var["gene_id"]]:
                db_var_match = db[var["gene_id"]]["premature_stop"]
            elif "large_deletion" in var["type"] and "large_deletion" in db[var["gene_id"]]:
                db_var_match = db[var["gene_id"]]["large_deletion"]
            if db_var_match:
                if "annotation" not in annotated_results["variants"][i]:
                    annotated_results["variants"][i]["annotation"] = {}
                for key in db_var_match:
                    annotated_results["variants"][i]["annotation"][key] = db_var_match[key]

    return annotated_results
