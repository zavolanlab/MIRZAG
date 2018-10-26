#!/bin/bash

mkdir -p logs/cluster_log
mkdir -p logs/local_log

snakemake \
	-s ../../../snakemake/Snakefile \
	-p \
	--configfile config.yml \
	--use-singularity \
	--singularity-args "--bind ${PWD},$PWD/../../,${PWD}/../../../" \
	--singularity-prefix ../../../../singularity/ \
	--cores 4
