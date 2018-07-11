##### BASE IMAGE #####
FROM python:2.7-slim

##### METADATA #####
LABEL base.image="python:2.7-slim"
LABEL version="1"
LABEL software="MIRZAG_scripts"
LABEL software.version="1"
LABEL software.description="Custom scripts for MIRZAG pipeline."
LABEL software.website="https://github.com/zavolanlab/MIRZAG"
LABEL software.documentation="https://github.com/zavolanlab/MIRZAG"
LABEL software.license="SPDX:GPL-3.0"
LABEL software.tags="Genomics, Transcriptomics"
LABEL maintainer="jaroslaw.surkont@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="SPDX:GPL-3.0"

##### INSTALL #####
COPY MIRZA/ /tmp/mirza/

RUN apt-get update && \
  apt-get install --yes \
    make \
    g++ \
    wget && \
  cd /tmp/mirza/ && \
  make && \
  mv MIRZA /usr/bin/ && \
  rm -rf /tmp/mirza && \
  cd /tmp && \
  wget http://contra.stanford.edu/contrafold/contrafold_v2_02.tar.gz && \
  tar xzf contrafold_*.tar.gz && \
  cd /tmp/contrafold/src && \
  sed -i '/^CXXFLAGS = / s/$/ -fpermissive/' Makefile && \
  sed -i '/#include <vector>/a #include <limits.h>' Utilities.hpp && \
  make && \
  mv contrafold /usr/bin/ && \
  rm -rf /tmp/contrafold* && \
  apt-get remove --purge --yes \
    make \
    g++ \
    wget && \
  apt-get autoremove --purge --yes

COPY requirements.txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt
COPY scripts/*.py /usr/bin/
