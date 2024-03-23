from pathogenprofiler.models import BarcodeResult, Variant, BamQC, FastaQC, DrVariant, GenomePosition
from .models import Lineage, TbDrVariant, TbVariant, ProfileResult, Spoligotype, LineageResult, Pipeline
from typing import List, Tuple , Union, Optional
from .utils import get_gene2drugs
import argparse
from pathogenprofiler.utils import shared_dict

def get_main_lineage(lineages: List[Lineage],max_node_skip: int = 1) -> Tuple[str, str]:
    """
    Get the main lineage and sublineage from a list of Lineage objects
    
    Arguments
    ---------
    lineages : List[Lineage]
        List of Lineage objects

    max_node_skip : int
        Maximum number of nodes to skip when collapsing lineages

    Returns
    -------
    Tuple[str, str]
        Tuple of main lineage and sublineage

    Examples
    --------
    >>> from tbprofiler import get_main_lineage
    >>> lineage = [
    ...     Lineage(fraction=0.8, lineage="lineage1", family="family1", spoligotype="spol1", rd="rd1"),
    ...     Lineage(fraction=0.2, lineage="lineage2", family="family2", spoligotype="spol2", rd="rd2"),
    ...     Lineage(fraction=0.2, lineage="lineage2.2", family="family2.2", spoligotype="spol2.2", rd="rd2.2"),
    ...     Lineage(fraction=0.8, lineage="lineage1.1", family="family1.1", spoligotype="spol1.1", rd="rd1.1"),
    ... ]
    >>> get_main_lineage(lineage)
    ('lineage1;lineage2', 'lineage1.1;lineage2.2')
    """
    def collapse_paths(paths):
        filtered_paths = []
        for p in sorted(paths,reverse=True):
            path_stored = any([p in x for x in filtered_paths])
            if not path_stored:
                filtered_paths.append(p)
        return filtered_paths

    def derive_path(x):
        return [".".join(x.split(".")[:i])for i in range(1,len(x.split(".")))] + [x]

    lin_freqs = {}
    pool = []
    for l in lineages:
        pool.append(l.lineage.replace("M.","M_"))
        lin_freqs[l.lineage.replace("M.","M_")] = float(l.fraction)
    routes = [";".join(derive_path(x)) for x in pool]
    paths = collapse_paths(routes)
    path_mean_freq = {}
    for path in paths:
        nodes = tuple(path.split(";"))
        nodes_skipped = sum([n not in pool for n in nodes])
        if nodes_skipped>max_node_skip: continue
        freqs = [lin_freqs[n] for n in nodes if n in lin_freqs]
        path_mean_freq[nodes] = sum(freqs)/len(freqs)
    main_lin = ";".join(sorted(list(set([x[0] for x in path_mean_freq])))).replace("_",".")
    sublin = ";".join(sorted(list(set([x[-1] for x in path_mean_freq])))).replace("_",".")
    return (main_lin,sublin)

def barcode2lineage(barcode: List[BarcodeResult]) -> List[Lineage]:
    """
    Convert a list of BarcodeResult objects to a list of Lineage objects
    
    Arguments
    ---------
    barcode : List[BarcodeResult]
        List of BarcodeResult objects
    
    Returns
    -------
    List[Lineage]
        List of Lineage objects
    
    Examples
    --------
    >>> from tbprofiler import barcode2lineage
    >>> barcode = [
    ...     BarcodeResult(id="barcode1", frequency=0.8, info=["family1", "spol1", "rd1"]),
    ...     BarcodeResult(id="barcode2", frequency=0.2, info=["family2", "spol2", "rd2"]),
    ... ]
    >>> barcode2lineage(barcode)
    [Lineage(fraction=0.8, lineage='barcode1', family='family1', spoligotype='spol1', rd='rd1'), Lineage(fraction=0.2, lineage='barcode2', family='family2', spoligotype='spol2', rd='rd2')]
    """

    lineage = []
    for d in barcode:
        lineage.append(
            Lineage(
                fraction=d.frequency,
                lineage=d.id,
                family=d.info[0],
                spoligotype=d.info[1],
                rd=d.info[2],
                support=d.support
            )
        )
    return lineage





def get_drtypes(dr_variants: List[TbDrVariant]) -> str:
    resistant_drugs = set()
    for var in dr_variants:
        resistant_drugs.update(var.get_drugs())

    FLQ_set = set(["levofloxacin","moxifloxacin","ciprofloxacin","ofloxacin"])
    groupA_set = set(["bedaquiline","linezolid"])

    rif = "rifampicin" in resistant_drugs
    inh = "isoniazid" in resistant_drugs
    flq = len(FLQ_set.intersection(resistant_drugs)) > 0
    gpa = len(groupA_set.intersection(resistant_drugs)) > 0

    if len(resistant_drugs)==0:
        drtype = "Sensitive"
    elif (rif and not inh) and not flq:
        drtype = "RR-TB"
    elif (inh and not rif):
        drtype = "HR-TB"
    elif (rif and inh) and not flq:
        drtype = "MDR-TB"
    elif rif and (flq and not gpa):
        drtype = "Pre-XDR-TB"
    elif rif and (flq and gpa):
        drtype = "XDR-TB"
    else:
        drtype = "Other"

    return drtype

unlist = lambda t: [item for sublist in t for item in sublist]


def variant_present(var,results):
    result = None
    if var['type']=='resistance_variant':
        for v in results['dr_variants']:
            if v['gene']==var['gene'] or v['locus_tag']==var['gene']:
                result = v
    else:
        for v in results['dr_variants'] + results['other_variants']:
            if (v['gene']==var['gene'] or v['locus_tag']==var['gene']) and v['type']==var['type']:
                result = v
    return result

def process_variants(
    variants: List[Union[Variant,DrVariant]], 
    bed_file: str
) -> List[Union[TbDrVariant,TbVariant]]:
    new_objects = []
    gene2drugs = get_gene2drugs(bed_file)
    for var in variants:
        dump = var.model_dump()
        dump['locus_tag'] = dump['gene_id']
        dump['gene_associated_drugs'] = gene2drugs.get(var.gene_name, [])
        if isinstance(var, DrVariant):
            new_objects.append(TbDrVariant(**dump))
        else:
            new_objects.append(TbVariant(**dump))
    return new_objects
        
def split_variants(
    variants: List[Variant],
    bed_file: str
) -> Tuple[List[TbDrVariant], List[TbVariant]]:
    variants = process_variants(variants, bed_file)

    dr_variants = []
    other_variants = []
    fail_variants = []
    for var in variants:
        if var.filter.upper() == "PASS":
            if isinstance(var, TbDrVariant):
                dr_variants.append(var)
            else:
                other_variants.append(var)
        else:
            fail_variants.append(var)
    return dr_variants,other_variants,fail_variants

def filter_missing_positions(missing_positions: List[GenomePosition]) -> List[GenomePosition]:
    for pos in missing_positions:
        who_annotations = [
            ann for ann in pos.annotation if 
                ann['type']=='who_confidence' and
                ann['confidence'] in ('Assoc w R - Interim','Assoc w R')
        ]
        other_annotations = [ann for ann in pos.annotation if ann['type']!='who_confidence']
        pos.annotation = who_annotations + other_annotations
    
    return [pos for pos in missing_positions if len(pos.annotation)>0]

def create_lineage_result(
    args: argparse.Namespace,
    lineage: List[Lineage]
):
    main_lineage, sub_lineage = get_main_lineage(lineage)
    data = {
        'id':args.prefix,
        'lineage':lineage,
        'sub_lineage':sub_lineage,
        'main_lineage':main_lineage,
        'tbprofiler_version':args.version,
        'db_version':args.conf['version'],
    }
    return LineageResult(**data)

def create_resistance_result(
    args: argparse.Namespace,
    notes: List[str],
    lineage: List[Lineage],
    spoligotype: Optional[Spoligotype],
    variants: List[Variant],
    qc: Union[BamQC, FastaQC]
) -> ProfileResult:
    dr_variants, other_variants, fail_variants = split_variants(variants,args.conf['bed'])
    main_lineage, sub_lineage = get_main_lineage(lineage)
    drtype = get_drtypes(dr_variants)
    pipeline = Pipeline(
        software_version=args.version,
        db_version=args.conf['version'],
        software=[{'process':k,'software':v} for k,v in shared_dict.items()]
    )
    if hasattr(qc, 'missing_positions'):
        qc.missing_positions = filter_missing_positions(qc.missing_positions)

    data = {
        'id':args.prefix,
        'notes':notes,
        'lineage':lineage,
        'spoligotype':spoligotype,
        'drtype':drtype,
        'dr_variants':dr_variants,
        'other_variants':other_variants,
        'qc_fail_variants':fail_variants,
        'sub_lineage':sub_lineage,
        'main_lineage':main_lineage,
        'pipeline':pipeline
    }

    return ProfileResult(**data, qc=qc)


def clean_up_duplicate_annotations(variants: Variant) -> None:
    """
    Remove duplicate annotations from a list of Variant objects
    
    Arguments
    ---------
    variants : List[Variant]
        List of Variant objects
    
    Returns
    -------
    None
    """
    confidence_levels = ['Assoc w R','Assoc w R - Interim','Uncertain significance','Not assoc w R - Interim','Not assoc w R']
    for var in variants:
        keys = set([(ann['type'],ann['drug']) for ann in var.annotation])
        new_annotations = []
        for key in keys:
            confidence_anns = []
            other_anns = []
            for ann in var.annotation:
                if ann['type']==key[0] and ann['drug']==key[1]:
                    if 'confidence' in ann and ann['confidence'] in confidence_levels:
                        confidence_anns.append(ann)
                    else:
                        other_anns.append(ann)
            # confidence_anns = [ann for ann in var.annotation if ann['type']==key[0] and ann['drug']==key[1] and 'confidence' in ann and ann['confidence'] in confidence_levels]
            # other_anns = [ann for ann in var.annotation if ann['type']==key[0] and ann['drug']==key[1] and ('confidence' not in ann or ann['confidence'] not in confidence_levels)]
            if len(confidence_anns)>0:
                confidence_anns = sorted(confidence_anns,key=lambda x: confidence_levels.index(x['confidence']))
                new_annotations.append(confidence_anns[0])

            new_annotations += other_anns
        var.annotation = new_annotations