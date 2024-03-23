# Webserver

## LSHTM TB-Profiler webserver

If you don't have access to a linux or macOS environment you can still use tb-profiler by using our webserver at http://tbdr.lshtm.ac.uk/

You can upload your next generation sequencing data in *fastQ* format. You can upload one or two (forward and reverse) fastq files. When you upload your data, the run will be be assigned a unique ID. Please take a note of this ID as you will need to to find your results later. Batch upload of samples is also possible.

## Setting up your own webserver

The code for the webserver is available [here](https://github.com/jodyphelan/tb-profiler-webserver). This is still in early development but you can use it to set up your own instance of the server. To do this run the following code:

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