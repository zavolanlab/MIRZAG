snakemake \
	--dag \
	-np \
	-s ../../snakemake/Snakefile \
	--configfile config.yml | dot -Tpng > dag.png
