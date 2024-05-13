import json
import os
from collections import defaultdict
from tqdm import tqdm
import logging
import csv
import argparse
from .models import ProfileResult
from typing import List, Tuple, Optional
from itol_config.interfaces import ColourStripConfigWriter, BinaryDataConfigWriter
from pydantic import BaseModel

class VariantDB:
    def __init__(self, json_db: Optional[dict] = None):
        self.samples2variants = defaultdict(set)
        self.variant2samples = defaultdict(set)
        self.variant_frequencies = {}
        self.samples = list()
        self.variant_rows = []
        if json_db:
            for gene in json_db:
                for mutation in json_db[gene]:
                        self.variant2samples[(gene,mutation)] = set()

    def add_result(self, result: ProfileResult) -> None:
        self.samples.append(result.id)
        for var in result.dr_variants + result.other_variants:
            key = (result.id,var.gene_name,var.change)
            self.variant_frequencies[key] = var.freq
            key = (var.gene_name,var.change)
            self.variant2samples[key].add(result.id)
            self.samples2variants[result.id].add(key)
            d = var.model_dump()
            d['sample'] = result.id
            d['drugs'] = ";".join([x['drug'] for x in var.drugs]) if hasattr(var,'drugs')>0 else "-"
            self.variant_rows.append(d)
    def get_frequency(self,key: Tuple[str,str,str]) -> float:
        return self.variant_frequencies.get(key,0.0)
    def get_variant_list(self) -> List[Tuple[str,str]]:
        return list(self.variant2samples.keys())
    def write_dump(self,filename: str) -> None:
        with open(filename,"w") as O:
            fields = ["sample","gene_name","change","freq","type","drugs"]
            writer = csv.DictWriter(O,fieldnames=fields)
            writer.writeheader()
            for row in self.variant_rows:
                d = {k:row[k] for k in fields}
                writer.writerow(d)

class TransmissionEdge(BaseModel):
    source: str
    target: str
    distance: float

    def dump(self):
        return {"source":self.source,"target":self.target,"properties":{"distance":self.distance}}
    
    def __eq__(self, __value: object) -> bool:
        if isinstance(__value, TransmissionEdge):
            if set([self.source, self.target]) == set([__value.source, __value.target]) and self.distance == __value.distance:
                return True
        return False
    def __hash__(self):
        sorted_samples = sorted([self.source,self.target])
        return hash((type(self),sorted_samples[0],sorted_samples[1],self.distance) )

def collate_results(args: argparse.Namespace) -> None:
    """
    Collate results from tbprofiler
    
    Arguments
    ---------
    args : argparse.Namespace
        Arguments from argparse
    """
    for d in args.dir:
        if not os.path.isdir(d):
            logging.error("\nERROR: Can't find directory %s\n" % d )
            exit()
    
    drugs = args.conf['drugs']

    samples = {}
    for d in args.dir:
        for f in os.listdir(f"{d}/"):
            if f.endswith(".results.json"):
                s = f.replace(".results.json","")
                samples[s] = f'{d}/{s}.results.json'
    
    if args.samples:
        samples = {s:samples[s] for s in [l.strip() for l in open(args.samples)]}

    rows = []
    edges = []
    variant_db = VariantDB(args.conf['json_db'])
    for s in tqdm(samples):
        data = json.load(open(samples[s]))
        res = ProfileResult(**data)
        variant_db.add_result(res)

        row = {
            "sample":s,
            "main_lineage":res.main_lineage,
            "sub_lineage":res.sub_lineage,
            "spoligotype":res.spoligotype.octal if res.spoligotype else "-",
            "drtype":res.drtype,
            "target_median_depth":res.qc.get_target_median_depth(),
            "pct_reads_mapped":res.qc.get_percent_reads_mapped(),
            "num_reads_mapped":res.qc.get_reads_mapped(),
            "num_dr_variants":len(res.dr_variants),
            "num_other_variants":len(res.other_variants),
        }
        for d in drugs:
            variants = [var.get_str() for var in res.dr_variants if d in var.get_drugs()]
            row[d] = ", ".join(sorted(variants)) if len(variants)>0 else "-"

        if args.mark_missing:
            drug_missing_variants = set()
            for var in res.qc.missing_positions:
                for ann in var.annotation:
                    if ann["type"]=="drug_resistance":
                        drug_missing_variants.add(ann['drug'])
            for d in drug_missing_variants:
                row[d] = "*%s" % row[d]

        rows.append(row)

        for linked_sample in res.linked_samples:
            edges.append(
                TransmissionEdge(
                    source=s,
                    target=linked_sample.sample,
                    distance=linked_sample.distance
                )
            )


    if args.format=="txt":
        args.sep = "\t"
    else:
        args.sep = ","

    outfile = f'{args.prefix}.{args.format}'
    with open(outfile,"w") as O:
        writer = csv.DictWriter(O,fieldnames=list(rows[0]),delimiter=args.sep)
        writer.writeheader()
        writer.writerows(rows)

    if args.itol:
        generate_itol_config(rows,drugs,args.prefix)

    generate_transmission_network(rows,edges,args.prefix)
    generate_distance_matrix(rows,edges,args.prefix)
    generate_variant_matrix(variant_db,args.prefix)

    variant_db.write_dump(f'{args.prefix}.variants.csv')


def generate_itol_config(rows: List[dict], drugs: list, prefix: str) -> None:
    """
    Generate itol configuration files

    Arguments
    ---------
    rows : List[dict]
        List of dictionaries containing results
    drugs : list
        List of drugs
    prefix : str
        Prefix for output files
    """
    all_lineage_cols = {"lineage1":"#104577","lineage2":"#ab2323","lineage3":"#18a68c","lineage4":"#f68e51","lineage5":"#7cb5d2","lineage6":"#fde05e","lineage7":"#bc94b7","lineage8":"#ccc9e7","lineage9":"#bd9391","Animal strains":"#f8e0c8","Other":"#000000","Not called": "#ffffff"}
    lineage_aggregation = {"": "Not called","M.caprae":"Animal strains","M.bovis":"Animal strains","M.orygis":"Animal strains"}
    lineage_dict = {r['sample']:lineage_aggregation.get(r["main_lineage"],r["main_lineage"]) if ";" not in r["main_lineage"] else "Other" for r in rows}
    lineages_present = set(lineage_dict.values())
    lineage_cols = {key:val for key,val in all_lineage_cols.items() if key in lineages_present}
    lineage_outfile = f'{prefix}.lineage.itol.txt'
    writer = ColourStripConfigWriter(lineage_dict,'Main Lineage',lineage_cols)
    writer.write(lineage_outfile)

    all_dr_cols = {"Sensitive":"#28a745","RR-TB":"#007bff","HR-TB":"#E0ACD5","MDR-TB":"#ffc107","Pre-XDR-TB":"#dc3545","XDR-TB":"#343a40","Other":"#f8f9fa"}
    dr_data = {row['sample']:row['drtype'] for row in rows}
    drtypes_present = set(dr_data.values())
    dr_cols = {key:val for key,val in all_dr_cols.items() if key in drtypes_present}
    drtypeoutfile = f'{prefix}.dr.itol.txt'
    writer = ColourStripConfigWriter(dr_data,'Drug Resistance',dr_cols)
    writer.write(drtypeoutfile)

    drug_binary_data = {row['sample']:{drug:1 if row[drug]!="-" else 0 for drug in drugs} for row in rows}
    drug_binary_outfile = f'{prefix}.dr.indiv.itol.txt'
    writer = BinaryDataConfigWriter(drug_binary_data,'Drug Resistance')
    writer.write(drug_binary_outfile)

def generate_transmission_network(rows: List[dict], edges: List[TransmissionEdge], prefix: str) -> None:
    """
    Generate transmission network in tjson format
    
    Arguments
    ---------
    rows : List[dict]
        List of dictionaries containing results
    edges : List[TransmissionEdge]
        List of edges in transmission network
    prefix : str
        Prefix for output files
    """
    row_lookup = {r['sample']:r for r in rows}
    if len(edges)>0:
        tmp_nodes = set()
        json_edges = []

        for e in edges:
            if e.source in row_lookup and e.target in row_lookup:
                tmp_nodes.add(e.source)
                tmp_nodes.add(e.target)
                json_edges.append(e.dump())
        nodes = []
        for n in tmp_nodes:
            nodes.append(
                {
                    "id":n,
                    "properties": {
                        "drtype":row_lookup[n]["drtype"],
                        "lineage":row_lookup[n]["main_lineage"],
                        "median_depth":row_lookup[n]["target_median_depth"]
                    }
                }
            )
        json.dump({"nodes":nodes,"edges":json_edges},open(prefix+".transmission_graph.json","w"))

def generate_distance_matrix(rows: List[dict], edges: List[TransmissionEdge], prefix: str) -> None:
    """
    Generate distance matrix for transmission network
    
    Arguments
    ---------
    rows : List[dict]
        List of dictionaries containing results
    edges : List[TransmissionEdge]
        List of edges in transmission network
    prefix : str
        Prefix for output files
    """
    samples = [r['sample'] for r in rows]
    if len(edges)>0:
        with open(prefix+".distance_matrix.txt","w") as MAT:
            MAT.write("samples\t%s\n" % "\t".join(samples))
            transformed_edges = {}
            for e in edges:
                transformed_edges[(e.source,e.target)] = e.distance
            for si in samples:
                row = [si]
                for sj in samples:
                    ss = tuple(sorted([si,sj]))
                    if si==sj:
                        row.append("0.0")
                    else:
                        row.append(str(transformed_edges[ss]) if ss in transformed_edges else "NA")
                MAT.write("%s\n" % "\t".join(row))


def generate_variant_matrix(variant_db: VariantDB, prefix: str) -> None:
    variant_list = variant_db.get_variant_list()
    with open(prefix+".variants.txt","w") as VAR:
        VAR.write("sample\t%s\n" % ("\t".join(["%s_%s" % (g,c) for g,c in variant_list])))
        for s in variant_db.samples:
            VAR.write("%s\t%s\n" % (s,"\t".join(["%.3f" % variant_db.get_frequency((s,gene,change)) for gene,change in variant_list])))

