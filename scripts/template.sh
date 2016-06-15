# Copyright (C) 2014 Rafal Gumienny <r.gumienny@unibas.ch>
#
# Distributed under terms of the GPL license.
#

# end when error
set -e
# raise error when variable is unset
set -u
# raise error when in pipe
set -o pipefail

# load specific modules
{%- if modules %}
{% for mymodule in modules -%}
module load {{ mymodule }}
{% endfor -%}
{% endif %}
# go to $TMPDIR
{% if copy -%}
cd {{ copydir }}

# Copy files to tmp directory
{% for from, to in copy.iteritems() -%}
cp -v {{ from }} {{ to }}
{% endfor %}
{% endif %}
#command to execute
{{ command }}

# move files back to server (usually output file)
{% if moveback -%}
{% for from, to in moveback.iteritems() -%}
mv {{ from }} {{ to }}
{% endfor %}
{% endif %}
