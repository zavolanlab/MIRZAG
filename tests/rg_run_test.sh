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

cwd=$(pwd)

# clean up
rm -rf config.ini
python ../MIRZA_G_pipeline.py clean -y

# substitute and run
sed "s~path_to_tests~${cwd}~g" config_template.ini > config.ini
python ../MIRZA_G_pipeline.py run --config config.ini
