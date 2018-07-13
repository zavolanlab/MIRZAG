#!/bin/bash

SINGULARITY_NOHTTPS=true snakemake -s ../../snakemake/Snakefile -p \
  --configfile config.yml --use-singularity --cores 4
