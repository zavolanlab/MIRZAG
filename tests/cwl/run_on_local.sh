#!/usr/bin/env bash

# SINGULARITY_NOHTTPS=true SINGULARITY_PULLFOLDER=${PWD} cwltool --singularity ../../cwl/mirzag_scan.cwl job.yml
SINGULARITY_NOHTTPS=true SINGULARITY_PULLFOLDER=${PWD} cwltool --singularity ../../cwl/mirzag_seed.cwl job.yml
