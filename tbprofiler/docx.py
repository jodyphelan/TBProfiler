import sys
from docxtpl import DocxTemplate
from collections import defaultdict
from .models import ProfileResult
from docx import Document
from typing import List
from copy import deepcopy
import logging
from abc import abstractmethod

class DocxResultTemplate:
    @abstractmethod
    def write_output(self, result: ProfileResult, conf: dict, outfile: str):
        """Write the output to the Word document"""
        pass
    
    @abstractmethod
    def get_template():
        """Return the path to the customised Word document template"""
        pass

def sanitize(d):
    d = d.replace("-","_")
    return d

def cache_cells(tab) -> List[List[str]]:
    _cells = tab._cells
    numcols = len(tab.columns)
    cells = []
    for row in tab.rows:
        rowcells = []
        for i in range(numcols):
            rowcells.append(_cells.pop(0))
        cells.append(rowcells)
    return cells

def merge_cells(filename: str) -> None:
    doc = Document(filename)
    
    def _merge_cells(tab, rows: List[int], column: int):
        
        if column >= len(tab.columns) - 1:
            return
        if len(rows)==1:
            return 
        c1 = tab.cell(rows[0], column)
        c2 = tab.cell(rows[-1], column)
        
        for r in rows[1:]:
            tab.cell(r, column).text = ""
        cm = c1.merge(c2)

        values_in_next_column = set([tab.rows[r].cells[column+1].text for r in rows])
        for val in values_in_next_column:
            rows_with_val = [r for r in rows if tab.c[r][column+1].text == val]
            _merge_cells(tab, rows_with_val, column+1)
    
    for tab in doc.tables:
        tab.c = cache_cells(tab)
        values_in_next_column = set([tab.rows[r].cells[0].text for r in range(1, len(tab.rows))])
        for val in values_in_next_column:
            rows_with_val = [r for r in range(1, len(tab.rows)) if tab.c[r][0].text == val]

            _merge_cells(tab, rows_with_val, 0)

    doc.save(filename)



def write_docx(result: ProfileResult,conf,outfile,template_file = None, plugin = None):
    if template_file is None:
        template_file = sys.prefix+"/share/tbprofiler/default_template.docx"
    
    if plugin:
        plugin_cls = plugin()
        variables = plugin_cls.write_output(result, conf,outfile)

    else:

        dr_variant_table = []
        comments = []
        
        for var in result.dr_variants:
            for drug in var.drugs:
                mutation = var.change if len(var.change) < 15 else var.change[:11]+"..."+var.change[-3:]
                if drug['comment'] not in comments and len(drug['comment']):
                    comments.append(drug['comment'])
                dr_variant_table.append({
                    'drug': drug['drug'],
                    'gene': var.gene_name,
                    'mutation': mutation,
                    'depth': var.depth,
                    'frequency': round(var.freq * 100,2),
                    'confidence': drug['confidence'],
                    'comment': comments.index(drug['comment'])+1 if len(drug['comment']) else ""
                })
        drugs_found = set([x['drug'] for x in dr_variant_table])
        for drug in conf['drugs']:
            if drug not in drugs_found:
                dr_variant_table.append({
                    'drug': drug,
                    'gene': "",
                    'mutation': "",
                    'depth': "",
                    'frequency': "",
                    'confidence': "",
                    'comment': ""
                })

        gene_qc = []
        if 'target_qc' in result.qc:
            for item in result.qc.target_qc:

                gene_qc.append({
                    'status': 'ok' if item.percent_depth_pass>0.9 else 'fail',
                    'gene': item.target,
                    'median_depth': item.median_depth,
                    'percent_depth_pass': item.percent_depth_pass,
                })

        dr_variant_table = sorted(dr_variant_table,key=lambda x: conf['drugs'].index(x['drug']))

        other_variants_table = []
        for var in result.other_variants:
            for ann in var.annotation:
                other_variants_table.append({
                    'gene': var.gene_name,
                    'mutation': var.change,
                    'depth': var.depth,
                    'frequency': round(var.freq * 100,2),
                    'drug': ann['drug'],
                    'confidence': ann['confidence'],
                })
        
        data = result.model_dump()



        variables = {
            'date':data['timestamp'].strftime("%d %b %Y"),
            'sublineage': data['sub_lineage'],
            'version': data['pipeline']['software_version'],
            'drtype': data['drtype'],
            'd': data,
            'dr_variants_table': dr_variant_table,
            'other_variants_table': other_variants_table,
            'comments': comments,
            'gene_qc': gene_qc
        }


        doc = DocxTemplate(template_file)
        doc.render(variables)
        doc.save(outfile)

        merge_cells(outfile)