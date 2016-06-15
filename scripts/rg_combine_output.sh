#! /bin/sh
#
# a.sh
# Copyright (C) 2014 Rafal Gumienny <r.gumienny@unibas.ch>
#
# Distributed under terms of the GPL license.
#
dir=$1

for i in $dir/*.tab
	do
		zcat $i
	done
