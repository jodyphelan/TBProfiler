# Windows

Sadly the bioinfomatic software that tb-profiler depends on does not run on Windows, so it can't be installed with conda. There is another option though - to run with docker. This quick guide will run through how to set up docker and set it up to run tb-profiler.

## Step 1 - Download and install docker

Go to https://docs.docker.com/desktop/windows/install/ to get the latest docker installation file. After it downloads, execute it and go through the recommended steps. During the installation process you might be asked to restart your computer a few times.

<img src="https://3546519222-files.gitbook.io/~/files/v0/b/gitbook-x-prod.appspot.com/o/spaces%2F-M9cvGy4eVqvGN5UqFAr%2Fuploads%2FxKJ8w7Z7ftFSPD8gcBHp%2FDocker_1.PNG?alt=media&token=c84c285b-e5df-4636-ba58-e6b1accb71d5">

After it finishes it will prompt you to install WSL 2.

<img src="https://3546519222-files.gitbook.io/~/files/v0/b/gitbook-x-prod.appspot.com/o/spaces%2F-M9cvGy4eVqvGN5UqFAr%2Fuploads%2Fn7JOnVUQzTIPQn50HXtQ%2FDocker_2.PNG?alt=media&token=9d734ef4-ec02-4ef2-a6bb-55056c597c72">

Again, follow the link on the prompt (https://aka.ms/wsl2kernal) and follow the instructions.

## Step 2 Install tb-profiler docker image

Each new release on conda is also built for docker and is available on quay.io. To get the image click the start menu and open up a program called PowerShell

<img src="https://3546519222-files.gitbook.io/~/files/v0/b/gitbook-x-prod.appspot.com/o/spaces%2F-M9cvGy4eVqvGN5UqFAr%2Fuploads%2F8wGPINCn7R80hxaG44Nn%2FDocker_8.jpg?alt=media&token=27877c76-fdf2-4718-a411-4b3646c5c5e0">

Once this is open type in the following:

```
docker pull quay.io/biocontainers/tb-profiler:4.3.0--pypyh5e36f6f_0
```

## Step 3 Create a container

Now that you have the image you can set up a container. To do this you can open up Docker Desktop which was installed as part of the Docker installation. Navigate to the "Images" tab and the hover over the tb-profiler image and click "Run". Expand the "Optional settings" section and fill in the Volumnes box. This allows you to map a folder on your local machine to a directory in the linux container:

* **Host path**: click the three dots to help you select the folder
* **Container path**: type in "/mnt"

<img src="https://3546519222-files.gitbook.io/~/files/v0/b/gitbook-x-prod.appspot.com/o/spaces%2F-M9cvGy4eVqvGN5UqFAr%2Fuploads%2FL5W92nGivltyk71fGRAY%2FDocker_5.PNG?alt=media&token=a30c07be-bdce-49b3-b137-df7e3c96c26b">

Then click "Run"

## Step 4 Open up a terminal from the container

Now navigate to the "Containers" tab in Docker Desktop, hover over the container and click on the "CLI" button as indicated below.

<img src="https://3546519222-files.gitbook.io/~/files/v0/b/gitbook-x-prod.appspot.com/o/spaces%2F-M9cvGy4eVqvGN5UqFAr%2Fuploads%2Fz8UyekpEUYAVAecsAZyy%2FDocker_6.PNG?alt=media&token=653ef919-ba37-4352-9c8a-0d0fb970442e">

This will open up a terminal window where you can run linux commands - including tb-profiler!

First navigate to the "/mnt" directory and then run tb-profiler as shown elsewhere in this documentation.

```
cd /mnt
tb-profiler profile -h
```