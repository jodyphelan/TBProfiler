# use the ubuntu base image
FROM ubuntu:14.04

RUN apt-get update && apt-get install -y \
	build-essential \
	autoconf \
	zlib1g-dev \
	libncurses-dev \
	libbz2-dev \
	libhdf5-dev \
	liblzma-dev \
	git \
	wget \
	unzip \
	python


RUN cd /opt \
	&& git clone --recursive https://github.com/jodyphelan/TBProfiler.git \
  && cd TBProfiler \
  && bash install_prerequisites.sh
