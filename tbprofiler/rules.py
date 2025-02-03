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

def inactivate_drug_resistance(variant: Variant):
    """
    Inactivate a drug resistance variant
    """
    for ann in variant.consequences[0].annotation:
        if ann['type']=='drug_resistance':
            ann['type'] = 'inactivated_drug_resistance'
    for csq in variant.consequences:
        for ann in csq.annotation:
            if ann['type']=='drug_resistance':
                ann['type'] = 'inactivated_drug_resistance'
            


class MmpR5WHORule(ProfilePlugin):
    """
    Epistasis rule for mmpL5/mmpR5
    """

    def process_variants(self,args,variants: List[Variant]):
        """Generic variant processing method"""

        v1 = search_variant(variants,drug='bedaquiline',gene_name='mmpR5')
        v2 = search_variant(variants,gene_name='mmpL5',type='LoF')

        v1_total_freq = math.ceil(sum([x.freq*100 for x in v1]))
        v1_total_freq = min(v1_total_freq,100)

        v2_total_freq = math.ceil(sum([x.freq*100 for x in v2]))
        v2_total_freq = min(v2_total_freq,100)

        if v1 and v2:
            v1_changes = ", ".join([x.change for x in v1])
            v2_changes = ", ".join([x.change for x in v2])
            note = f"Loss of function mutation(s) detected in mmpL5 ({v2_changes}) which may abrogate the effect of the genetically linked mmpR5 mutation(s) ({v1_changes})."
            
            if v2_total_freq==100:
                inactivate_drug_resistance(v1)
            freq_diff = v1_total_freq - v2_total_freq
            if freq_diff > 10:
                note += f" However, the combined frequency of the mmpR5 mutation(s) is {freq_diff}% higher than the mmpL5 mutation(s), indicating a potential resistant subpopulation. Please consult the raw data for more information."
            
            args.notes.append(note)


class eisWHORule(ProfilePlugin):
    """
    Epistasis rule for mmpL5/mmpR5
    """

    def process_variants(self,args,variants: List[Variant]):
        """Generic variant processing method"""

        v1 = search_variant(variants,drug=['kanamycin','amikacin'],gene_name='eis')
        v2 = search_variant(variants,gene_name='eis',type='LoF')

        v1_total_freq = math.ceil(sum([x.freq*100 for x in v1]))
        v1_total_freq = min(v1_total_freq,100)

        v2_total_freq = math.ceil(sum([x.freq*100 for x in v2]))
        v2_total_freq = min(v2_total_freq,100)

        if v1 and v2:
            v1_changes = ", ".join([x.change for x in v1])
            v2_changes = ", ".join([x.change for x in v2])
            note = f"Loss of function mutation(s) detected in eis ({v2_changes}) which may abrogate the effect of the genetically linked eis promoter mutation(s) ({v1_changes})."

            if v2_total_freq==100:
                inactivate_drug_resistance(v1)
            freq_diff = v1_total_freq - v2_total_freq
            if freq_diff > 10:
                note += f" However, the combined frequency of the eis promoter mutation(s) is {freq_diff}% higher than the eis coding mutation(s), indicating a potential resistant subpopulation. Please consult the raw data for more information."
            
            args.notes.append(note)


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

