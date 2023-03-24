# Fichiers de sortie

## Formats de sortie

### Sortie par défaut 

Par défaut, un fichier de sortie au format `.json` est produit dans un répertoire appelé `results`. Ce fichier contient des informations telles que les mutations trouvées, la lignée ainsi que quelques métriques de contrôle de qualité (QC). Ce format est parfait pour être chargé dans un script pour un traitement en aval, mais il n'est pas très lisible par l'homme.

### Sorties de texte

Pour un format plus lisible par l'homme, vous pouvez utiliser les options `--csv` et `--txt` pour générer des sorties au format csv et texte également. Ces fichiers contiendront plusieurs tableaux listant des informations similaires au format json. Pour les utilisateurs avancés, il est possible de personnaliser ce format en fournissant un fichier modèle. Ce modèle doit être écrit dans le [jinja templating language](https://jinja.palletsprojects.com/en/3.1.x/). Une variable nommée `d` sera disponible dans le modèle. Cette variable contient exactement la même structure que le fichier json. Par exemple, un modèle simple pourrait être :

```
Custom TB-Profiler report

Lineage: {{d['sublin']}}
Drug-resistance: {{d['drtype']}}
```

### Sortie Docx

Il est également possible de produire un beau rapport au format docx qui peut être visualisé dans Word ou converti en pdf. L'avantage de ce format est qu'il peut contenir des images, du texte formaté, etc. Pour créer ce rapport, un fichier modèle doit être fourni. Comme pour le format texte personnalisé, le modèle docx doit contenir des variables jinja qui seront remplies avec les données de l'échantillon lorsqu'un rapport est généré. Des détails sur les variables disponibles pour le moteur de création de modèles peuvent être trouvés [ici](https://github.com/jodyphelan/TBProfiler/blob/master/tbprofiler/docx.py).

![example_report](https://github.com/jodyphelan/TBProfiler/raw/docs/docs/assets/img/report_example.png)

## Générer de fichiers de synthèse

Les résultats d'un grand nombre d'exécutions peuvent être rassemblés dans un seul tableau à l'aide de la commande suivante :

```
tb-profiler collate
```

Cela créera automatiquement un certain nombre de fichiers de résultats groupés à partir de tous les fichiers de résultats individuels dans le répertoire résultats. Si vous souhaitez générer ce fichier pour un sous-ensemble des runs, vous pouvez fournir une liste avec les noms des runs en utilisant l'option `--samples`. Le préfixe des fichiers de sortie est tbprofiler par défaut, mais il peut être modifié avec l'option `--prefix`.

### Rédiger vos propres scripts de résumé

La fonction collate extrait les mutations de résistance aux médicaments et la lignée, mais il se peut que vous souhaitiez extraire d'autres caractéristiques présentes dans les fichiers de résultats json individuels. J'ai créé un petit tutoriel sur la façon de procéder [ici](https://jodyphelan.gitbook.io/tb-profiler/writing-a-custom-collate-script).

### Fichiers iTOL

Plusieurs fichiers sont produits par la fonction `tb-profile collate`. Parmi ceux-ci se trouvent plusieurs fichiers de configuration qui peuvent être utilisés avec [iTOL](http://itol.embl.de/) pour annoter les arbres phylogénétiques. Un petit arbre et des fichiers de configuration ont été placés dans le répertoire example_data. Pour l'utiliser, naviguez sur le site web d'iTOL et téléchargez le fichier tbprofiler.tree en utilisant le bouton upload sur la barre de navigation. Une fois le fichier téléchargé, vous accéderez à une visualisation de l'arbre. Pour ajouter l'annotation, glissez-déposez le fichier sur l'arbre dans le navigateur. Vous devriez maintenant voir une figure similaire à celle ci-dessous. Les annotations suivantes sont incluses :

* Lignée
* Classes de résistance aux médicaments (sensible, résistant aux médicaments, MDR, XDR)
* Appels de résistance aux médicaments pour les médicaments individuels, les cercles remplis représentent la résistance.

<img href="https://github.com/jodyphelan/TBProfiler/raw/docs/docs/docs/assets/images/itol_example.png">
