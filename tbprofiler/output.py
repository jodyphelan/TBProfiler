import json
import time
def load_drgus(dr_bed_file):
    drugs = set()
    for l in open(dr_bed_file):
        arr = l.rstrip().split()
        for d in arr[5].split(";"):
            drugs.add(d)
    return sorted(list(drugs))


def write_results_2(self):
    drugs = load_drgus(self.params["dr_bed_file"])
    o = open(self.params["txt_results"],"w")
    o.write("#### %(prefix)s results ####\n\n" % self.params)
    o.write("#### Drug resistance associated small variants ####\n")
    for d in drugs:
        dr_vars = [x for x in self.small_dr_variants["dr"] if x["drug"]==d]
        res = "R" if len(dr_vars)>0 else ""
        muts = ["%(gene)s (%(change)s, %(freq)s)" % x for x in dr_vars]
        o.write("%s\t%s\t%s\n" % (d,res,", ".join(muts)))
    dr_drugs = [x["drug"] for x in self.small_dr_variants["dr"]]

    MDR = "R" if ("ISONIAZID" in dr_drugs and "RIFAMPICIN" in dr_drugs) else ""
    XDR = "R" if MDR=="R" and ( "AMIKACIN" in dr_drugs or "KANAMYCIN" in dr_drugs or "CAPREOMYCIN" in dr_drugs ) and ( "FLUOROQUINOLONES" in dr_drugs) else ""
    o.write("MDR\t%s\n" % MDR)
    o.write("XDR\t%s\n" % XDR)

    o.write("\n\n#### Big Deletions in candidate genes ####\n")
    for l in self.deletions:
        o.write("\t".join([str(l[x]) for x in ["drug","gene","len","good_cov",'float_gene_present']])+"\n")

    o.write("\n\n#### Lineage ####\n")
    o.write("\n".join(["%s\t%s\t%s\t%s\t%s" % (x["lin"],x["frac"],x["family"],x["spoligotype"],x["rd"]) for x in self.lineage]))

    o.write("\n\n#### Other small variants in candidate genes ####\n")
    for l in self.small_dr_variants["other"]:
        o.write("\t".join([str(l[x]) for x in ["chr","genome_pos","gene","locus_tag","change","type","freq"]])+"\n")
    o.close()

    json_dict = {"small_variants_other":self.small_dr_variants["other"],"small_variants_dr":self.small_dr_variants["dr"],"del":self.deletions,"lineage":self.lineage,"id":self.params["prefix"]}

    json.dump(json_dict,open(self.params["json_results"],"w"))

def write_results_1(self):
    drugs = load_drgus(self.params["dr_bed_file"])
    o = open(self.params["txt_results"],"w")
    o.write("#### %(prefix)s results ####\n\n" % self.params)
    o.write("#### Drug resistance associated small variants ####\n")
    for l in self.small_dr_variants["dr"]:
        o.write("\t".join([str(l[x]) for x in ["drug","chr","genome_pos","gene","locus_tag","change","type","freq"]])+"\n")
    o.write("\n\n#### Big Deletions in candidate genes ####\n")
    for l in self.deletions:
        o.write("\t".join([str(l[x]) for x in ["drug","gene","len","good_cov",'float_gene_present']])+"\n")
    o.write("\n\n#### Lineage ####\n")
    o.write("\n".join(["%s\t%s\t%s\t%s\t%s" % (x["lin"],x["frac"],x["family"],x["spoligotype"],x["rd"]) for x in self.lineage]))
    o.write("\n\n#### Other small variants in candidate genes ####\n")
    for l in self.small_dr_variants["other"]:
        o.write("\t".join([str(l[x]) for x in ["chr","genome_pos","gene","locus_tag","change","type","freq"]])+"\n")
    o.close()

    json_dict = {"small_variants_other":self.small_dr_variants["other"],"small_variants_dr":self.small_dr_variants["dr"],"del":self.deletions,"lineage":self.lineage,"id":self.params["prefix"]}

    json.dump(json_dict,open(self.params["json_results"],"w"))

def write_tex(self):
    self.drugs = set()
    for l in open(self.params["dr_bed_file"]):
        arr = l.rstrip().split()
        for d in arr[5].split(";"):
            self.drugs.add(d)


    write_json(self)

    json_results = json.load(open(self.params["json_results"]))
    results = {}
    drug_table = []
    for x in json_results["small_variants_dr"]:
        for d in x["drug"].split(";"):
            if d not in results: results[d] = set()
            results[d].add("\\textit{%s} %s" % (x["gene"],x["change"]))
    for d in self.drugs:
        if d in results:
            results[d] = ", ".join(results[d]) if len(results[d])>0 else ""
            r = "R" if len(results[d])>0 else ""
        else:
            results[d] = ""
            r = ""
        drug_table.append({"Drug":d,"Genotypic Resistance":r,"Mutations":results[d]})

    pipeline_tbl = [{"Analysis":"Mapping","Program":self.params["mapper"] if self.params["mapping"] else "N/A"},{"Analysis":"Variant Calling","Program":self.caller}]

    o = open(self.params["tex_results"],"w")
    tex_strings = {}
    tex_strings["id"] = self.params["prefix"]
    tex_strings["date"] = time.ctime()
    tex_strings["strain"] = json_results["sublin"]
    tex_strings["drtype"] = json_results["drtype"]
    tex_strings["dr_report"] = dict_list2tex(drug_table,["Drug","Genotypic Resistance","Mutations"])
    tex_strings["lineage_report"] = dict_list2tex(json_results["lineage"],["lin","frac","family","spoligotype","rd"],{"lin":"Lineage","frac":"Estimated fraction"})
    tex_strings["other_var_report"] = dict_list2tex(json_results["small_variants_other"],["genome_pos","locus_tag","change","freq"],{"genome_pos":"Genome Position","locus_tag":"Locus Tag","freq":"Estimated fraction"})
    tex_strings["pipeline"] = dict_list2tex(pipeline_tbl,["Analysis","Program"])
    tex_strings["version"] = self.tbpv
    o.write(load_tex(tex_strings))
    o.close()






def write_json(self):
    json_dict = {"small_variants_other":self.small_dr_variants["other"],"small_variants_dr":self.small_dr_variants["dr"],"del":self.deletions,"lineage":self.lineage,"id":self.params["prefix"]}
    json_dict["main_lin"] = sorted([x["lin"] for x in json_dict["lineage"]])[0] if len(json_dict["lineage"])>0 else "-"
    json_dict["sublin"] = sorted([x["lin"] for x in json_dict["lineage"]])[-1] if len(json_dict["lineage"])>0 else "-"
    dr_drugs = [x["drug"] for x in json_dict["small_variants_dr"]]
    self.dr_drugs = dr_drugs
    MDR = "R" if ("ISONIAZID" in dr_drugs and "RIFAMPICIN" in dr_drugs) else "-"
    XDR = "R" if MDR=="R" and ( "AMIKACIN" in dr_drugs or "KANAMYCIN" in dr_drugs or "CAPREOMYCIN" in dr_drugs ) and ( "FLUOROQUINOLONES" in dr_drugs) else "-"
    drtype = "Sensitive"
    if XDR=="R":
        drtype="XDR"
    elif MDR=="R":
        drtype="MDR"
    elif len(dr_drugs)>0:
        drtype="Drug-resistant"
    json_dict["XDR"] = XDR
    json_dict["MDR"] = MDR
    json_dict["drtype"] = drtype
    json.dump(json_dict,open(self.params["json_results"],"w"))

def load_tex(tex_strings):
    return r"""
%%%% LyX 2.3.0 created this file.  For more info, see http://www.lyx.org/.
\documentclass[english]{report}
\renewcommand{\familydefault}{\sfdefault}
\usepackage[T1]{fontenc}
\usepackage[latin9]{inputenc}
\usepackage[a4paper]{geometry}
\geometry{verbose,tmargin=2cm,bmargin=2cm,lmargin=2cm,rmargin=2cm}
\setcounter{secnumdepth}{3}
\setcounter{tocdepth}{3}
\usepackage{graphicx}

\makeatletter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LyX specific LaTeX commands.
%% Because html converters don't know tabularnewline
\providecommand{\tabularnewline}{\\}

\makeatother

\usepackage{babel}
\begin{document}

\chapter*{TBProfiler report}

The following report has been generated by TBProfiler.
\section*{Summary}
\begin{description}
\item [{ID:}] %(id)s
\item [{Date:}] %(date)s
\item [{Strain:}] %(strain)s
\item [{Drug-resitsance:}] %(drtype)s
\end{description}

\section*{Lineage report}

%(lineage_report)s

\section*{Resistance report}

%(dr_report)s

\section*{Other variants report}

%(other_var_report)s

\section*{Analysis pipeline specifications}

\begin{description}
\item [{Version:}] %(version)s
\end{description}

%(pipeline)s

\section*{Disclaimer}

This tool is for \textbf{Research Use Only} and is offered free for
use. The London School of Hygiene and Tropical Medicine shall have
no liability for any loss or damages of any kind, however sustained
relating to the use of this tool.

\section*{Citation}

Coll, F. \textit{et al}. Rapid determination of anti-tuberculosis
drug resistance from whole-genome sequences. \textit{Genome Medicine}
7, 51. 2015

\end{document} """ % tex_strings


def dict_list2tex(l,columns = None, mappings = {}):
    headings = l[0].keys() if not columns else columns
    rows = []
    header = " & ".join([mappings[x].title() if x in mappings else x.title() for x in headings])+"\\tabularnewline"
    for row in l:
        r = " & ".join(["%.3f" % row[x] if type(row[x])==float else str(row[x]).replace("_", " ") for x in headings])
        rows.append(r)
    column_def = "".join(["l" for _ in headings])
    str_rows = "\\tabularnewline\n".join(rows)+"\\tabularnewline"
    return "\\begin{tabular}{%s}\n%s\n\\hline\n%s\n\\hline\n\end{tabular}" % (column_def,header,str_rows)
