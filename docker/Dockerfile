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
COPY requirements.txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt
COPY scripts/*.py /usr/bin/
