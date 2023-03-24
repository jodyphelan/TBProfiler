# Serveur web

## LSHTM TB-Profiler serveur web

Si vous n'avez pas accès à un environnement linux ou macOS, vous pouvez toujours utiliser tb-profiler en utilisant notre serveur web à l'adresse http://tbdr.lshtm.ac.uk/.

Vous pouvez télécharger vos données de séquençage de nouvelle génération au format *fastQ*. Vous pouvez télécharger un ou deux fichiers fastq (forward et reverse). Lorsque vous téléchargez vos données, un identifiant unique est attribué à la série. Veuillez noter cet identifiant, car vous en aurez besoin pour retrouver vos résultats ultérieurement. Il est également possible de télécharger des échantillons par lots.

## Mise en place de votre propre serveur web

Le code du serveur web est disponible [ici](https://github.com/jodyphelan/tb-profiler-webserver). Il s'agit encore d'un développement précoce, mais vous pouvez l'utiliser pour mettre en place votre propre instance du serveur. Pour ce faire, exécutez le code suivant :

```
# Install libraries
git clone https://github.com/jodyphelan/tb-profiler-webserver.git
cd tb-profiler-webserver
python setup.py install

# Run flask
export FLASK_APP=tbprofiler_web
export FLASK_ENV=development
flask run

# Run rabbit-mq server
rabbitmq-server

# Run celery
celery -A tbprofiler_web.worker worker --loglevel=info --concurrency=1
```
