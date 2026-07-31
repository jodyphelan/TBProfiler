from .plugins import ProfilePlugin
import logging
from pathogenprofiler.models import Variant
from .models import ProfileResult
from typing import List
import math
import argparse
import re

SILENT_MUTATIONS = ('synonymous_mutation','initiator_codon_variant','stop_retained_variant','start_retained_variant')

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

class lineageResistanceRule(Rule):
    """
    A rule which adds intrinsic resistance to report
    if a lineage is found
    """
    __domain__ = 'result'
    def process_result (
            self,
            args: argparse.Namespace,
            lineage: str,
            drug: str,
            result: ProfileResult,
            **kwargs
        ):
        for l in result.lineage:
            if l.lineage==lineage:
                result.notes.append(kwargs['note'])

        
class SequencingArtefectRule(Rule):
    """
    A rule which flags a warding when a variant is found in a certain
    position in a gene
    """
    __domain__ = 'result'
    def process_result (
        self,
        args: argparse.Namespace,
        result: ProfileResult,
        **kwargs
    ):
        for var in result.dr_variants + result.other_variants:
            if var.gene_name==kwargs['gene'] or var.gene_id==kwargs['gene']:
                r = re.search('c.([0-9]+).>.',var.nucleotide_change)
                if r:
                    if r.group(1) == str(kwargs['gene_position']):
                        logging.debug(f"Adding sequencing artefact note for {var.gene_name} {var.change} at position {r.group(1)}")

                        result.notes.append(kwargs['note'])
    
        

class epistasisRule(Rule):
    """
    Epistasis rule
    """
    __domain__ = 'variants'
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
        annotation_ids = set()
        d = args.conf['json_db']
        for g in d:
            for m in d[g]:
                for a in d[g][m]['annotations']:
                    for key in a:
                        annotation_ids.add(key)

        for var in result.other_variants:
            confidence = {}
            for ann in var.annotation:
                if 'confidence' in ann:
                    confidence[ann['drug']] = ann['confidence']
            
            for drug in var.gene_associated_drugs:
                if drug not in confidence:
                    if var.type in SILENT_MUTATIONS:                        
                        confidence[drug] = 'Not Assoc W R - Interim'
                    else:
                        confidence[drug] = 'Uncertain significance'
                    ann = {
                        'type':'who_confidence',
                        'drug':drug,
                        'confidence':confidence[drug],
                        'comment':'Not found in WHO catalogue'
                    }
                    for key in annotation_ids:
                        if key not in ann:
                            ann[key] = ''
                    var.annotation.append(ann)
                    logging.debug(f'{var.gene_name} {var.change} does not have a confidence value for {drug}. Setting it to {confidence[drug]}')

class CrossResistanceRule(Rule):

    """
    Cross resisrance rule rule
    """
    __domain__ = 'variants'
    def process_variants(
            self, 
            source_drug: str,
            target_drug: str,
            note: str,
            args: argparse.Namespace,
            variants: List[Variant], 
            **kwargs
        ):
        for var in variants:
            for ann in list(var.annotation):
                new_ann = ann.copy()
                if 'drug' in ann and source_drug == ann['drug']:
                    new_ann['drug'] = target_drug
                    new_ann['comment'] = note
                    var.annotation.append(new_ann)
                    logging.debug(f'Adding cross resistance annotation to {var.gene_name} {var.change} for {target_drug} based on {source_drug} annotation')

def is_resistance_variant(var: Variant, drug: str) -> bool:
    for ann in var.annotation:
        if ann['type']=='drug_resistance' and ann['drug']==drug:
            return True
    return False

class ResistanceLevelRule(Rule):
    """
    Docstring for ResistanceLevelRule
    """
    __domain__ = 'variants'
    def process_variants(
            self, 
            drug: str,
            target_variants: List[dict],
            args: argparse.Namespace,
            variants: List[Variant], 
            **kwargs
        ):

        high_level_resistance_variants = []
        for v in target_variants:
            if v['resistance_level'].lower() == 'high':
                high_level_resistance_variants.append((v['gene_name'],v['change']))
        for var in variants:
            key = (var.gene_name, var.change)
            if is_resistance_variant(var,drug):
                    resistance_level = 'high' if key in high_level_resistance_variants else kwargs['default_resistance_level']
                    var.annotation.append({
                        'type':'resistance_level',
                        'drug':drug,
                        'resistance_level': resistance_level,
                        'comment':'High level resistance mutation' if key in high_level_resistance_variants else 'Low level resistance mutation'
                    })
                    logging.debug(f'Setting resistance level to {resistance_level} for {var.gene_name} {var.change} for {drug}')

class CompensatoryRule(Rule):
    """
    Docstring for CompensatoryRule
    """
    __domain__ = 'variants'
    def process_variants(
        self, 
        args: argparse.Namespace,
        variants: List[Variant], 
        resistance_gene: str,
        compensatory_gene: str,
        compensatory_mutations: List[str],
        drug: str,
        original_confidence: str,
        updated_confidence: str,
        note: str,
        **kwargs
    ):
        compensatory_variant_present = False
        for var in variants:
            if var.gene_name==compensatory_gene and var.change in compensatory_mutations:
                compensatory_variant_present = True
                logging.debug(f"Found compensatory mutation {var.gene_name} {var.change} which may abrogate the effect of a linked resistance mutation in {resistance_gene} for {drug}")
                break
        
        high_level_resistance_variants = False
        for var in variants:
            resistance_variant = False
            for ann in var.annotation:
                if ann['type']=='drug_resistance' and ann['drug']==drug and ann['confidence']=="Assoc w R":
                    resistance_variant = True
                    break

            if var.gene_name==resistance_gene and resistance_variant:
                for ann in var.annotation:
                    if ann['type']=='who_confidence' and ann['drug']==drug and ann['confidence']==original_confidence:
                        high_level_resistance_variants = True
                        logging.debug(f"Found resistance mutation {var.gene_name} {var.change} with confidence {original_confidence} for {drug} which may be abrogated by compensatory mutation(s) in {compensatory_gene}")
                        break

        logging.debug(f"Compensatory variant present: {compensatory_variant_present}"
                      f"\nHigh level resistance variants present: {high_level_resistance_variants}")
        if compensatory_variant_present and not high_level_resistance_variants:
            for var in variants:
                if var.gene_name==resistance_gene:
                    for ann in var.annotation:
                        if ann['type']=='who_confidence' and ann['drug']==drug and ann['confidence']==original_confidence:
                            logging.debug(f"Updating confidence from {original_confidence} to {updated_confidence} for {var.gene_name} {var.change} for {drug} based on presence of compensatory mutation(s) in {compensatory_gene}")
                            ann['confidence'] = updated_confidence
                            ann['type'] = 'drug_resistance'
                            ann['comment'] = note
    
