# Mise à jour

TB-Profiler fait l'objet d'un développement rapide et constant. Si vous envisagez d'utiliser le programme dans le cadre de votre travail, assurez-vous d'utiliser la version la plus récente ! De même, la base de données n'est pas statique et fait l'objet d'améliorations constantes ; assurez-vous donc d'utiliser la version la plus récente. Si vous utilisez TBProfiler dans votre travail, veuillez indiquer la version de l'outil et de la base de données, car ils sont développés indépendamment l'un de l'autre.

### Mise à jour de la base de données

De nouvelles mutations/gènes sont régulièrement ajoutés à la base de données. Exécutez les opérations suivantes pour vous assurer que vous êtes à jour.

```
tb-profiler update_tbdb
```

### Re-profilage rapide

Si vous disposez d'une nouvelle base de données de mutations mais qu'aucun nouveau gène n'a été ajouté, vous pouvez rapidement reprofiler vos échantillons en exécutant la procédure suivante.

```
tb-profiler reprofile /path/to/result.json
```

Cela peut être utile lorsque vous avez ajouté quelques mutations vous-même ou que vous êtes sûr qu'aucun nouveau gène n'a été ajouté dans la mise à jour. Si vous n'êtes pas sûr, il est plus prudent d'exécuter à nouveau l'étape de profilage complet.

