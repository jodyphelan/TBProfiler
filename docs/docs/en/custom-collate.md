# Custom collate scripts

Result files for each individual run are by default stored in json files (ending in **.results.json**). These can be collated into a merged file to display the drug resistance mutations using `tb-profiler collate`. However, you may want to customise what information is extracted from the result files and displayed in the collated output. This page will provide an outline how to write your own parser.

## Script skeleton

Open up your favourite editor, paste the following code and save the file as `tbprofiler_custom_script.py`

```
#! /usr/bin/env python

import sys
import csv
import json
import argparse
import os
from collections import defaultdict
from tqdm import tqdm
import pathogenprofiler as pp
import tbprofiler as tbprofiler

def main(args):
    bed_file = "%s/share/tbprofiler/%s.bed" % (sys.base_prefix,args.db)
    locus_tag2drugs = tbprofiler.get_lt2drugs(bed_file)

    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]

    for s in tqdm(samples):
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))


parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--dir',default="results/",type=str,help='Directory containing results')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
```

This will form the basic structure of the parser. Let's see what each section does.

**Lines 3-11**: Here we import libraries that are useful in the script. 

**Lines 13**: This line starts the definition of the mainfunction. It takes one argument in the form of an object of arguments which is constructed using the argument parser.

**Line 14**: The location of the bed file is stored in the variable bed_file. This is constructed using the prefix of the database supplied in the command line arguments (tbdb by default).

**Line 15**: We use the `tbprofiler.get_lt2drugs` function to extract the locus_tag to associated drug mapping in a dictionary with the following format:

```
{
    "Rv0667":["rifampicin"],
    "Rv1484":["isoniazid","ethionamide"],
}
```

**Line 17-20**: If the `--samples` argument was supplied, the file is read in and stored as a list. Otherwise the list is derived from looking the the result directory (specified with `--dir`).

**Lines 22-23**: A loop is constructed to loop over the samples. For each iteration the json result file is loaded into a dictionary structure. The key-value structure is explained below.

**Lines 26-34**: These lines parse the arguments supplied to the script (for example `--samples`). You can add in new arguments and they will appear as properties of the args object in the main function.

### Result data structure

The individual result files are read into a dictionary on line 23. There is a lot of information in the result. Below is a snippet of the dictionary.

```
{
    "main_lin": "lineage4",
    "sublin": "lineage4.9",
    "dr_variants": [
        {
            "chrom": "Chromosome",
            "genome_pos": 761155,
            "ref": "C",
            "alt": "T",
            "freq": 1,
            "feature_id": "CCP43410",
            "type": "missense_variant",
            "nucleotide_change": "c.1349C>T",
            "protein_change": "p.Ser450Leu",
            "alternate_consequences": [],
            "change": "p.Ser450Leu",
            "locus_tag": "Rv0667",
            "gene": "rpoB",
            "drugs": [
                {
                    "type": "drug",
                    "drug": "rifampicin",
                    "literature": "10.1128/AAC.01093-18",
                    "confers": "resistance"
                }
            ],
            "annotation": [
                {
                    "type": "resistance_association_confidence",
                    "drug": "rifampicin",
                    "confidence": "high"
                }
            ]
        }
    ],
    "other_variants": [
        {
            "chrom": "Chromosome",
            "genome_pos": 7362,
            "ref": "G",
            "alt": "C",
            "freq": 1,
            "feature_id": "CCP42728",
            "type": "missense_variant",
            "nucleotide_change": "c.61G>C",
            "protein_change": "p.Glu21Gln",
            "alternate_consequences": [],
            "change": "p.Glu21Gln",
            "locus_tag": "Rv0006",
            "gene": "gyrA",
            "gene_associated_drugs": [
                "ofloxacin",
                "levofloxacin",
                "fluoroquinolones",
                "ciprofloxacin",
                "moxifloxacin"
            ]
        },
    ]
}
```

The mutations are split into two lists "dr_variants" and "other_variants". The dr_variants are all associated with drug resistance. The other_variants are mutations in drug resistant genes but have not (yet!) been associated with resistance. Some of the other mutations may be phylogenetic mutations which have nothing to do with resistance while others may be but might not have enough evidence yet.

At the moment if you run the script and point it to a directory with some result files in using `--dir, it will read in all the data but not actually do anything with it. 

## Working example script

As an example I will extend its functionality to output a simple CSV file with the name of the sample and the drug-resistance type.

First I will create a CSV writer object. Add the following code before the for loop.

```
OUT = open(args.out,"w")
writer = csv.DictWriter(OUT, fieldnames = ["sample","dr-class"])
writer.writeheader()
```

Close the file handle after the loop finishes using

```
OUT.close()
```

We have used `args.out` property when creating the filehandle. For this to work we have to add an extra argument to the argument parser. Add the following for code together with the other arguments:

```
parser.add_argument('--out', type=str, help='Name of CSV output', required=True)
```

Now lets add in the code that determines the drug resistance class. By default tb-profiler only classes samples into sensitive, MDR, XDR and drug-resistant. Let's extend this to find pre-MDR and pre-XDR.

First we can create a set which we will then populate with a list of the drugs which the sample is resistant to. We can do this by looping over the variants. Add in the following code inside of the for loop:

```
resistant_drugs = set()
for var in data["dr_variants"]:
    for d in var["drugs"]:
        resistant_drugs.add(d["drug"])
```

Drug class definitions

Samples can be classed into different types using the following definitions

| Type      | Drugs resistance |
|-----------|------------------|
| Sensitive | No drug resistance |
| Pre-MDR          |   Rifampicin **or** isonisazid               |
|  MDR | Rifampicin **and** isoniazid |
| Pre-XDR | MDR **and** any fluoriquinolone | 
| XDR | MDR **and** (any fluoriquinolone **and** any group A drug) |
| Other | Resistance to any drug but none of the above categories |

### Coding this up

We can see straight away that the definitions are just a few boolean operations. We can store the boolean values by looking up the drugs in the `resistant_drugs` set we created before. First we will create a list of the fluoroquinolones (FLQ) and second line injectables (SLI). Add this code anywhere in the main function but before the loop:

```
FLQ_set = set(["moxifloxacin","levofloxacin","ciprofloxacin","ofloxacin"])
GPA_set = set(["bedaquiline","linezolid"])
```

Then create the boolean values. Put these inside the loop after you have created the data variable:

```
resistant_drugs = set()
for var in data["dr_variants"]:
    for d in var["drugs"]:
        resistant_drugs.add(d["drug"])
rif = "rifampicin" in resistant_drugs
inh = "isoniazid" in resistant_drugs
flq = len(FLQ_set.intersection(resistant_drugs)) > 0
gpa = len(GPA_set.intersection(resistant_drugs)) > 0
```

To test for rifampicin and isoniazid we simple check if they exists in resistant_drugs. To check for FLQ and GPA resistance we find the intersections between `FLQ_set` and resistant_drugs and count. If the value is greater than 0 then we have resistance in that category.

Now we are ready to perform the construction of the classes using our boolean variables. We can run through the table above using a series of  if-else statements.

```
if len(resistant_drugs)==0:
    drtype = "Sensitive"
elif (rif and not inh) or (inh and not rif):
    drtype = "Pre-MDR"
elif (rif and inh) and (not flq and not sli):
    drtype = "MDR"
elif (rif and inh) and ( flq and not gpa ):
    drtype = "Pre-XDR"
elif (rif and inh) and (flq and gpa):
    drtype = "XDR"
else:
    drtype = "Other"
```

We have stored the drug resistance class in the variable `drtype`. Now all that is left to do is to write this to the CSV file. We can do this using:

```
writer.writerow({"sample":s, "dr-class":drtype})
```

Your script should now look like this: 

```
#! /usr/bin/env python

import sys
import csv
import json
import argparse
import os
from collections import defaultdict
from tqdm import tqdm
import pathogenprofiler as pp
import tbprofiler as tbprofiler

def main(args):
    bed_file = "%s/share/tbprofiler/%s.bed" % (sys.base_prefix,args.db)
    locus_tag2drugs = tbprofiler.get_lt2drugs(bed_file)

    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]

    FLQ_set = set(["moxifloxacin","levofloxacin","ciprofloxacin","ofloxacin"])
    GPA_set = set(["bedaquiline","linezolid"])

    OUT = open(args.out,"w")
    writer = csv.DictWriter(OUT, fieldnames = ["sample","dr-class"])
    writer.writeheader()

    for s in tqdm(samples):
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
        resistant_drugs = set()
        for var in data["dr_variants"]:
            for d in var["drugs"]:
                resistant_drugs.add(d["drug"])
        rif = "rifampicin" in resistant_drugs
        inh = "isoniazid" in resistant_drugs
        flq = len(FLQ_set.intersection(resistant_drugs)) > 0
        gpa = len(GPA_set.intersection(resistant_drugs)) > 0

        if len(resistant_drugs)==0:
            drtype = "Sensitive"
        elif (rif and not inh) or (inh and not rif):
            drtype = "Pre-MDR"
        elif (rif and inh) and (not flq and not gpa):
            drtype = "MDR"
        elif (rif and inh) and ( flq and not gpa ):
            drtype = "Pre-XDR"
        elif (rif and inh) and (flq and gpa):
            drtype = "XDR"
        else:
            drtype = "Other"

        writer.writerow({"sample":s, "dr-class":drtype})

    OUT.close()

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--out', type=str, help='Name of CSV output', required=True)
parser.add_argument('--samples', type=str, help='File with samples')
parser.add_argument('--dir', default="results/", type=str, help='Directory containing results')
parser.add_argument('--db', default="tbdb", type=str, help='Database name')
parser.add_argument('--suffix', default=".results.json", type=str, help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
```

Try running using it and check if the output looks like it should:

```
python tbprofiler_custom_script.py --dir /path/to/dir --out test.csv
```

If it is working as expected, try adding in other pieces of info to the output. For example, you can extract the sub-lineage by accessing `data["sublin"]`. Have a look at the result file structure ( you can open it with a web browser to render it in a more readable format). 

Hopefully you can use this template as a boilerplate for your own uses. There are a few examples available [here](https://github.com/jodyphelan/tbprofiler_scripts) if you are interested. If you have any questions or if you have a cool script that you would like to make available to the community don't hesitate to contact me by email or through [GitHub](https://github.com/jodyphelan/TBProfiler/issues). 
