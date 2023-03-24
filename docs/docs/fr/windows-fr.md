# Windows

Malheureusement, le logiciel bioinfomatic dont dépend tb-profiler ne fonctionne pas sous Windows, il ne peut donc pas être installé avec conda. Il y a cependant une autre option : l'exécuter avec docker. Ce guide rapide explique comment installer docker et le configurer pour exécuter tb-profiler.

## Étape 1 - Télécharger et installer docker

Allez sur https://docs.docker.com/desktop/windows/install/ pour obtenir le dernier fichier d'installation de Docker. Une fois le fichier téléchargé, exécutez-le et suivez les étapes recommandées. Au cours du processus d'installation, il se peut que l'on vous demande de redémarrer votre ordinateur à plusieurs reprises.

<img src="https://3546519222-files.gitbook.io/~/files/v0/b/gitbook-x-prod.appspot.com/o/spaces%2F-M9cvGy4eVqvGN5UqFAr%2Fuploads%2FxKJ8w7Z7ftFSPD8gcBHp%2FDocker_1.PNG?alt=media&token=c84c285b-e5df-4636-ba58-e6b1accb71d5">

Une fois l'opération terminée, vous serez invité à installer WSL 2.

<img src="https://3546519222-files.gitbook.io/~/files/v0/b/gitbook-x-prod.appspot.com/o/spaces%2F-M9cvGy4eVqvGN5UqFAr%2Fuploads%2Fn7JOnVUQzTIPQn50HXtQ%2FDocker_2.PNG?alt=media&token=9d734ef4-ec02-4ef2-a6bb-55056c597c72">

Là encore, suivez le lien sur l'invite (https://aka.ms/wsl2kernal) et suivez les instructions.

## Étape 2 Installer l'image docker de tb-profiler

Chaque nouvelle version de conda est également construite pour docker et est disponible sur quay.io. Pour obtenir l'image, cliquez sur le menu Démarrer et ouvrez un programme appelé PowerShell.

<img src="https://3546519222-files.gitbook.io/~/files/v0/b/gitbook-x-prod.appspot.com/o/spaces%2F-M9cvGy4eVqvGN5UqFAr%2Fuploads%2F8wGPINCn7R80hxaG44Nn%2FDocker_8.jpg?alt=media&token=27877c76-fdf2-4718-a411-4b3646c5c5e0">

Une fois ouvert, tapez ce qui suit :

```
docker pull quay.io/biocontainers/tb-profiler:4.3.0--pypyh5e36f6f_0
```

## Étape 3 Créer un container

Maintenant que vous avez l'image, vous pouvez configurer un container. Pour ce faire, vous pouvez ouvrir Docker Desktop qui a été installé dans le cadre de l'installation de Docker. Naviguez vers l'onglet "Images" et survolez l'image tb-profiler et cliquez sur "Run". Développez la section "Optional settings" et remplissez la case Volumnes. Cela vous permet de faire correspondre un dossier sur votre machine locale à un répertoire dans le container linux :

* **Host path**: cliquez sur les trois points pour vous aider à sélectionner le dossier
* **Container path**: ecrir "/mnt"

<img src="https://3546519222-files.gitbook.io/~/files/v0/b/gitbook-x-prod.appspot.com/o/spaces%2F-M9cvGy4eVqvGN5UqFAr%2Fuploads%2FL5W92nGivltyk71fGRAY%2FDocker_5.PNG?alt=media&token=a30c07be-bdce-49b3-b137-df7e3c96c26b">

Cliquez ensuite sur "Run"

## Étape 4 Ouvrir un terminal à partir du conteneur

Naviguez maintenant vers l'onglet "Containers" dans Docker Desktop, survolez le conteneur et cliquez sur le bouton "CLI" comme indiqué ci-dessous.

<img src="https://3546519222-files.gitbook.io/~/files/v0/b/gitbook-x-prod.appspot.com/o/spaces%2F-M9cvGy4eVqvGN5UqFAr%2Fuploads%2Fz8UyekpEUYAVAecsAZyy%2FDocker_6.PNG?alt=media&token=653ef919-ba37-4352-9c8a-0d0fb970442e">

Cela ouvrira une fenêtre de terminal où vous pourrez exécuter des commandes linux - y compris tb-profiler !

Naviguez d'abord jusqu'au répertoire "/mnt", puis exécutez tb-profiler comme indiqué ailleurs dans cette documentation.

```
cd /mnt
tb-profiler profile -h
```
