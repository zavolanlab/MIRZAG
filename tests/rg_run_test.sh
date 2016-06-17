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

sed "s~path_to_tests~${cwd}~g" config_template.ini
