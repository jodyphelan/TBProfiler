from .models import DrugResistanceInterpretation, IntrinsicDrugResistance, DrVariant
from typing import List
import json

def get_dr_interpretation(
    drugs: List[str],
    intrinsic_dr:List[IntrinsicDrugResistance], 
    dr_variants: List[DrVariant],
    other_variants: list
) -> List[DrugResistanceInterpretation]:
    interpretation_dict = {}
    for drug in drugs:
        interpretation_dict[drug] = {
            'drug': drug,
            'resistance_mechanisms': [],
            'resistance_level': [],
            'warnings': [],
            'comments': set()
        }
    for mechanism in intrinsic_dr:
        if mechanism.drug in interpretation_dict:
            interpretation_dict[mechanism.drug]['resistance_mechanisms'].append('intrinsic_resistance')
            interpretation_dict[mechanism.drug]['comments'].add(mechanism.comment)

    for mechanism in dr_variants:
        for d in mechanism.drugs:
            drug = d['drug']
            if drug in interpretation_dict:
                interpretation_dict[drug]['resistance_mechanisms'].append(f'{mechanism.gene_name}:{mechanism.change}')
                interpretation_dict[drug]['comments'].add(d['comment'])
                if 'WARNING' in d['comment']:
                    interpretation_dict[drug]['warnings'].append(d['comment'])
                if d['resistance_level']:
                    interpretation_dict[drug]['resistance_level'].append(d['resistance_level'].lower())

    for drug in interpretation_dict:
        interpretation_dict[drug]['comments'] = list(interpretation_dict[drug]['comments'])
    
    for drug in interpretation_dict:
        d = interpretation_dict[drug]
        interpretation = ""
        if len(d['warnings']) > 0 and len(d['resistance_mechanisms'])==len(d['warnings']):
            interpretation = "W"
        elif "high" in d['resistance_level']:
            interpretation = "R (high)"
        elif d['resistance_level'].count('low') >1:
            interpretation = "R (high)"
        elif d['resistance_level'].count('low')==1:
            interpretation = "R (low)"
        elif len(d['resistance_mechanisms']) > 0:
            interpretation = "R"

        d['interpretation'] = interpretation
    interpretation_list = [DrugResistanceInterpretation(**d) for d in interpretation_dict.values()]

    return interpretation_list
          
            

    
