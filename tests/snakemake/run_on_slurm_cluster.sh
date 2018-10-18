#!/bin/bash

mkdir -p logs/cluster_log
mkdir -p logs/local_log

snakemake \
	-s ../../snakemake/Snakefile \
	-p \
	--configfile config.yml \
	--cores 20 \
	--local-cores 2 \
	--cluster-config cluster.json \
	--cluster "sbatch --cpus-per-task={cluster.threads} --mem={cluster.mem} --qos={cluster.queue} --time={cluster.time} --job-name={cluster.name} -o {cluster.out} --workdir . -p scicore" \
	--use-singularity \
	--singularity-args "--bind ${PWD},$PWD/../,${PWD}/../../" \
	--singularity-prefix ../../../singularity/










