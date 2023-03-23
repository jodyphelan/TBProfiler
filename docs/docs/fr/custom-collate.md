# Scénarios d'assemblage personnalisés

Les fichiers de résultats pour chaque analyse individuelle sont par défaut stockés dans des fichiers json (terminant par **.results.json**). Ils peuvent être rassemblés dans un fichier fusionné pour afficher les mutations de résistance aux médicaments en utilisant `tb-profiler collate`. Cependant, il se peut que vous souhaitiez personnaliser les informations extraites des fichiers de résultats et affichées dans le fichier fusionné. Cette page fournit un aperçu de la façon d'écrire votre propre analyseur syntaxique.

## Squelette du script

Ouvrez votre éditeur préféré, collez le code suivant et enregistrez le fichier sous le nom `tbprofiler_custom_script.py`

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

Cela constituera la structure de base de l'analyseur. Voyons ce que fait chaque section.

**Lignes 3-11**: Ici, nous importons des bibliothèques qui seront utiles dans le script. 

**Lignes 13**: Cette ligne commence la définition de la fonction principale. Elle prend un argument sous la forme d'un objet d'arguments construit à l'aide de l'analyseur d'arguments.

**Ligne 14** : L'emplacement du fichier bed est stocké dans la variable bed_file. Celle-ci est construite en utilisant le préfixe de la base de données fournie dans les arguments de la ligne de commande (tbdb par défaut).

**Ligne 15** : Nous utilisons la fonction `tbprofiler.get_lt2drugs` pour extraire la correspondance entre le locus_tag et le médicament associé dans un dictionnaire au format suivant :

```
{
    "Rv0667":["rifampicin"],
    "Rv1484":["isoniazid","ethionamide"],
}
```

**Ligne 17-20** : Si l'argument `--samples` a été fourni, le fichier est lu et stocké sous forme de liste. Sinon, la liste est obtenue en regardant le répertoire de résultat (spécifié avec `--dir`).

**Lignes 22-23** : Une boucle est construite pour parcourir les échantillons. Pour chaque itération, le fichier résultats json est chargé dans une structure de dictionnaire. La structure d'une valeur clé est expliquée ci-dessous.

**Lignes 26-34** : Ces lignes analysent les arguments fournis au script (par exemple `--samples`). Vous pouvez ajouter de nouveaux arguments et ils apparaîtront en tant que propriétés de l'objet args dans la fonction principale.

### Structure des données de résultat

Les fichiers de résultats individuels sont lus dans un dictionnaire à la ligne 23. Le résultat contient beaucoup d'informations. Voici un extrait du dictionnaire ci-dessous.

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

Les mutations sont réparties en deux listes : "dr_variants" et "other_variants". Les dr_variants sont tous associés à la résistance aux médicaments. Les autres_variants sont des mutations dans des gènes de résistance aux médicaments, mais n'ont pas (encore !) été associées à la résistance. Certaines des autres mutations peuvent être des mutations phylogénétiques qui n'ont rien à voir avec la résistance, tandis que d'autres peuvent l'être mais n'ont pas encore été suffisamment prouvées.

Actuellement, si vous exécutez le script et que vous le dirigez vers un répertoire contenant des fichiers résultats en utilisant `--dir, il lira toutes les données mais n'en fera rien.

## Exemple de script opérationnel

À titre d'exemple, j'étendrai sa fonctionnalité pour produire un simple fichier CSV contenant le nom de l'échantillon et le type de résistance aux médicaments.

Tout d'abord, je vais créer un objet CSV writer. Ajoutez le code suivant avant la boucle for.

```
OUT = open(args.out,"w")
writer = csv.DictWriter(OUT, fieldnames = ["sample","dr-class"])
writer.writeheader()
```

Fermer la poignée du fichier après la fin de la boucle en utilisant


```
OUT.close()
```
Nous avons utilisé la propriété `args.out` lors de la création de la poignée de fichier. Pour que cela fonctionne, nous devons ajouter un argument supplémentaire à l'analyseur d'arguments. Ajoutez le code suivant avec les autres arguments :

```
parser.add_argument('--out', type=str, help='Name of CSV output', required=True)
```

Ajoutons maintenant le code qui détermine la classe de résistance aux médicaments. Par défaut, tb-profiler ne classe les échantillons que dans les catégories suivantes : sensible, MDR, XDR et résistant aux médicaments. Étendons cette classification pour trouver les classes pré-MDR et pré-XDR.

Nous pouvons tout d'abord créer un ensemble que nous alimenterons ensuite avec une liste des médicaments auxquels l'échantillon est résistant. Nous pouvons le faire en parcourant en boucle les variantes. Ajoutez le code suivant à l'intérieur de la boucle for :

```
resistant_drugs = set()
for var in data["dr_variants"]:
    for d in var["drugs"]:
        resistant_drugs.add(d["drug"])
```

Définitions des classes de médicaments

Les échantillons peuvent être classés en différents types à l'aide des définitions suivantes

| Type      | Drugs resistance |
|-----------|------------------|
| Sensitive | No drug resistance |
| Pre-MDR          |   Rifampicin **or** isonisazid               |
|  MDR | Rifampicin **and** isoniazid |
| Pre-XDR | MDR **and** any fluoriquinolone | 
| XDR | MDR **and** (any fluoriquinolone **and** any group A drug) |
| Other | Resistance to any drug but none of the above categories |

### Coder cela

Nous pouvons voir tout de suite que les définitions sont juste quelques opérations booléennes. Nous pouvons stocker les valeurs booléennes en recherchant les médicaments dans l'ensemble `resistant_drugs` que nous avons créé précédemment. Tout d'abord, nous allons créer une liste des fluoroquinolones (FLQ) et des médicaments injectables de seconde ligne (SLI). Ajoutez ce code n'importe où dans la fonction principale, mais avant la boucle :

```
FLQ_set = set(["moxifloxacin","levofloxacin","ciprofloxacin","ofloxacin"])
GPA_set = set(["bedaquiline","linezolid"])
```

Créez ensuite les valeurs booléennes. Placez-les dans la boucle après avoir créé la variable de données :

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

Pour tester la rifampicine et l'isoniazide, il suffit de vérifier s'ils existent dans resistant_drugs. Pour vérifier la résistance au FLQ et à la GPA, nous trouvons les intersections entre `FLQ_set` et resistant_drugs et nous les comptons. Si la valeur est supérieure à 0, nous avons une résistance dans cette catégorie.

Nous sommes maintenant prêts à construire les classes à l'aide de nos variables booléennes. Nous pouvons parcourir le tableau ci-dessus à l'aide d'une série d'instructions if-else.

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

Nous avons stocké la classe de résistance aux médicaments dans la variable `drtype`. Il ne reste plus qu'à l'écrire dans le fichier CSV. Nous pouvons le faire en utilisant :

```
writer.writerow({"sample":s, "dr-class":drtype})
```

Votre script devrait maintenant ressembler à ceci : 

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

Essayez de l'utiliser et vérifiez si le résultat est correct :

```
python tbprofiler_custom_script.py --dir /path/to/dir --out test.csv
```

Si cela fonctionne comme prévu, essayez d'ajouter d'autres informations à la sortie. Par exemple, vous pouvez extraire le sous-lignage en accédant à `data["sublin"]`. Jetez un oeil à la structure du fichier résultat (vous pouvez l'ouvrir avec un navigateur web pour le rendre dans un format plus lisible).

Nous espérons que vous pourrez utiliser ce modèle comme un modèle standard pour vos propres utilisations. Il y a quelques exemples disponibles [ici](https://github.com/jodyphelan/tbprofiler_scripts) si vous êtes intéressés. Si vous avez des questions ou si vous avez un script sympa que vous aimeriez mettre à la disposition de la communauté, n'hésitez pas à me contacter par email ou via [GitHub](https://github.com/jodyphelan/TBProfiler/issues). 

