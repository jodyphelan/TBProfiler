from .plugins import ProfilePlugin
import logging
from pathogenprofiler.models import Variant
from .models import ProfileResult
from typing import List
import math
import argparse


def search_variant(variants: List[Variant], **kwargs) -> List[Variant]:
    type_expansions = {
        'lof': ('frameshift_variant','stop_gained','transcript_ablation','feature_ablation')
    }
    found_variants = set()
    for var in variants:
        for csq in var.consequences:
            tests = {}
            for key,val in kwargs.items():

                if isinstance(val,str):
                    values = [val]
                else:
                    values = val

                for ele in values:
                    if ele.lower() in type_expansions:
                        values = values + list(type_expansions[ele.lower()])

                test = False
                if (hasattr(csq,key) and vars(csq)[key] in values):
                    test = True
                for ann in csq.annotation:
                    if key in ann and ann[key] in values:
                        test = True

                tests[key] = test

            if all(tests.values()):
                found_variants.add(var)

    return list(found_variants)

def inactivate_drug_resistance(variants: List[Variant]):
    """
    Inactivate a drug resistance variant
    """
    for var in variants:
        for csq in var.consequences:
            for ann in csq.annotation:
                if ann['type']=='drug_resistance':
                    ann['type'] = 'inactivated_drug_resistance'
            
class Rule(ProfilePlugin):
    pass

class epistasisRule(Rule):
    """
    Epistasis rule
    """

    def process_variants(
            self, 
            source:dict, 
            target:dict, 
            args: argparse.Namespace,
            variants: List[Variant], 
            source_inactivation_freq_cutoff:int=100, 
            target_escape_freq_cutoff:int=10,
            **kwargs
        ):
        """Generic variant processing method"""

        source_vars = search_variant(variants,**source)
        target_vars = search_variant(variants,**target)

        source_vars_total_freq = math.ceil(sum([x.freq*100 for x in source_vars]))
        source_vars_total_freq = min(source_vars_total_freq,100)

        target_vars_total_freq = math.ceil(sum([x.freq*100 for x in target_vars]))
        target_vars_total_freq = min(target_vars_total_freq,100)
        
        if source_vars and target_vars:
            source_vars_changes = ", ".join([x.change for x in source_vars])
            target_vars_changes = ", ".join([x.change for x in target_vars])
            note = f"Mutation(s) detected in {source_vars[0].gene_name} ({source_vars_changes}) which may abrogate the effect of the genetically linked {target_vars[0].gene_name} mutation(s) ({target_vars_changes})."
            
            if source_vars_total_freq>=source_inactivation_freq_cutoff:
                inactivate_drug_resistance(target_vars)
            freq_diff = target_vars_total_freq - source_vars_total_freq 
            if freq_diff > target_escape_freq_cutoff:
                note += f" However, the combined frequency of the {target_vars[0].gene_name} mutation(s) is {freq_diff}% higher than the {source_vars[0].gene_name} mutation(s), indicating a potential resistant subpopulation."
            
            note += " Please consult the raw data for more information."
            args.notes.append(note)

def apply_epistasis_rule(args: argparse.Namespace, variants: List[Variant], parameters: dict):
    logging.debug(f"Applying epistasis rule with parameters: {parameters}")
    epistasisRule().process_variants(args=args,variants=variants,**parameters)


class SetConfidence(ProfilePlugin):

    def process_result(self, args: argparse.Namespace, result: ProfileResult):
        for var in result.other_variants:
            confidence = {}
            for ann in var.annotation:
                if 'confidence' in ann:
                    confidence[ann['drug']] = ann['confidence']
            
            for drug in var.gene_associated_drugs:
                if drug not in confidence:
                    if var.type=='synonymous_mutation':
                        confidence[drug] = 'Not Assoc W R - Interim'
                    else:
                        confidence[drug] = 'Uncertiain Significance'
                    var.annotation.append(
                        {
                            'type':'who_confidence',
                            'drug':drug,
                            'confidence':confidence[drug],
                            'comment':''
                        }
                    )
                    logging.debug(f'{var.gene_name} {var.change} does not have a confidence value for {drug}. Setting it to {confidence[drug]}')

