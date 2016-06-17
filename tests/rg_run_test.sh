#! /bin/sh
#
# rg_run_test.sh
# Copyright (C) 2016 Rafal Gumienny <r.gumienny@unibas.ch>
#
# Distributed under terms of the GPL license.
#


# end when error
set -e
# raise error when variable is unset
set -u
# raise error when in pipe
set -o pipefail

if [ "$#" -lt 1  ]; then
	echo "Usage: rg_run_test.sh clean/run [MIRZA/binary/path] ['CONTRAfold/binary/path']"
	exit
fi


cwd=$(pwd)

if [ $1 == 'clean' ]
then
	echo "Cleaning"
	rm -rf config.ini
	python ../MIRZA_G_pipeline.py clean -y
elif [ $1 == 'run' ]
then
	echo "Running tests"
	if [ -z ${2+x} ]; then
		echo "MIRZA is not set: setting default (MIRZA)."
		mirza="MIRZA"
	else
		mirza=$2
	fi

	if [ -z ${3+x} ]; then
		echo "CONTRAfold is not set: setting default (contrafold)."
		contrafold="contrafold"
	else
		contrafold=$3
	fi
	rm -rf config.ini
	python ../MIRZA_G_pipeline.py clean -y
	sed "s~path_to_tests~${cwd}~g" config_template.ini > config.ini
	sed -i "s~MIRZApath~${mirza}~g" config.ini
	sed -i "s~contrafoldpath~${contrafold}~g" config.ini
	python ../MIRZA_G_pipeline.py run --config config.ini
elif [ $1 == 'help' ]
then
	echo "Usage: rg_run_test.sh clean/run [MIRZA/binary/path] ['CONTRAfold/binary/path']"
else
	echo "Usage: rg_run_test.sh clean/run [MIRZA/binary/path] ['CONTRAfold/binary/path']"
fi
