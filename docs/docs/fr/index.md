# Introduction

TB-Profiler est un outil en ligne de commande utilisé pour traiter les données de séquençage du génome entier afin de prédire la lignée et la résistance aux médicaments. Le pipeline recherche les petits variants et les grandes délétions associés à la résistance aux médicaments. Il indique également la lignée. Par défaut, il utilise trimmomatic pour découper les reads, BWA (ou minimap2 pour nanopore) pour aligner les sequences reads au génome de référence et bcftools (ou GATKv4/freebayes) pour identifier les variants. 

<img src="https://files.gitbook.com/v0/b/gitbook-legacy-files/o/assets%2F-M9cvGy4eVqvGN5UqFAr%2F-M9dZ5yUBJa3XVGD-wRl%2F-M9d_diTzl9Ae2kLM03R%2Ftb-profiler_uml.svg?alt=media&token=7cf9db30-a0c6-448b-a750-f4a041f34478">
