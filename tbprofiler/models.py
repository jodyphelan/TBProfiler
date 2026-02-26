from pydantic import BaseModel, Field
from typing import Optional, List, Union
from pathogenprofiler import object_list2text, dict_list2text
from pathogenprofiler.models import BamQC, FastaQC, VcfQC, Variant, DrVariant, BarcodePosition
from datetime import datetime

__model_schema_version__ = '1.0.0'


class Lineage(BaseModel):
    """
    A class to hold information about a lineage
    
    Attributes
    ----------
    fraction : float
        Fraction of reads belonging to this lineage
    lineage : str
        Name of the lineage
    family : str
        Family names associated with the lineage
    rd : Optional[str]
        RDs associated with the lineage
    """
    fraction: float
    lineage: str
    family: str
    rd: Optional[str] = None
    support: List[BarcodePosition]

class Pipeline(BaseModel):
    """
    A class to hold information about the TB-Profiler pipeline
    
    Attributes
    ----------
    tbprofiler_version : str
        TB-Profiler version
    db_version : dict
        TB-Profiler database version
    software : List[dict]
        Software used in the pipeline
    """
    software_version: str
    db_version: dict
    software: List[dict]

class Result(BaseModel):
    """
    A class to hold information about a TBProfiler result
    
    Attributes
    ----------
    id : str
        Sample ID
    timestamp : datetime
        Time of analysis
    tbprofiler_version : str
        TBProfiler version
    db_version : dict
        TBProfiler database version
    """
    schema_version: str = __model_schema_version__
    id: str
    timestamp: datetime = Field(default_factory=datetime.now)
    pipeline: Pipeline

class TbDrVariant(DrVariant):
    locus_tag: str
    gene_associated_drugs: List[str] = []

class TbVariant(Variant):
    locus_tag: str
    gene_associated_drugs: List[str] = []

class Spacer(BaseModel):
    name: str
    seq: str
    count: int

class Spoligotype(BaseModel):
    binary: str
    octal: str
    family: Optional[str]
    SIT: Optional[str]
    countries: Optional[str]
    spacers: List[Spacer]

    def __repr__(self) -> str:
        return self.octal

class LinkedSample(BaseModel):
    sample: str
    distance: float
    positions: List[int]

class ProfileResult(Result):
    notes: List[str] = []
    lineage: Optional[List[Lineage]] = []
    main_lineage: str = None
    sub_lineage: str = None
    spoligotype: Optional[Spoligotype] = None
    drtype: str
    dr_variants: List[TbDrVariant] = []
    other_variants: List[TbVariant] = []
    qc_fail_variants: List[Union[TbDrVariant,TbVariant]] = []
    qc: Union[BamQC, FastaQC, VcfQC]
    linked_samples: List[LinkedSample] = []
    gene_name2locus_tag: dict = {}

    def get_qc(self, sep="\t"):
        if isinstance(self.qc, (BamQC, FastaQC)):
            lt2gene = {v: k for k, v in self.gene_name2locus_tag.items()}
            rows = []
            for item in self.qc.target_qc:
                if hasattr(item, "model_dump"):
                    row = item.model_dump()
                elif isinstance(item, dict):
                    row = dict(item)
                else:
                    row = vars(item).copy()

                target = row.get("target", "")
                if target in lt2gene:
                    row["locus_tag"] = target
                    row["gene_name"] = lt2gene[target]
                else:
                    row["locus_tag"] = self.gene_name2locus_tag.get(target, row.get("locus_tag", ""))
                    row["gene_name"] = target

                ordered_row = {
                    "locus_tag": row.get("locus_tag", ""),
                    "gene_name": row.get("gene_name", row.get("target", "")),
                }
                for key, val in row.items():
                    if key in ("locus_tag", "gene_name", "target"):
                        continue
                    ordered_row[key] = val
                rows.append(ordered_row)
            text = dict_list2text(l=rows, sep=sep)
        else:
            text = "Not available for VCF input"
        return text

    def get_missing_pos(self,sep="\t"):
        if isinstance(self.qc, (BamQC,)):
            text = object_list2text(
                self.qc.missing_positions,
                mappings={
                    "pos":"Genome Position",
                    "annotation.locus_tag":"Locus Tag",
                    "annotation.gene_name":"Gene name",
                    "annotation.variant":"Variant",
                    "annotation.drug":"Drug",
                    "annotation.source":"Source",
                    "annotation.confidence":"Confidence",
                    "depth":"Depth"
                },
                sep=sep
            )
        else:
            text = "Not available for input data type"
        return text

class LineageResult(Result):
    lineage: Optional[List[Lineage]] = []
    main_lineage: str = None
    sub_lineage: str = None

    def get_lineage(self):
        if self.lineage:
            return object_list2text(l = self.lineage)
        else:
            return "Not available"
