#!/bin/bash
set -e

jobcard="{{ jobcard }}"
export infile="{{ infile }}"
export ngenerate="{{ ngenerate }}"

{% for folder in replica_folders %}
export OUTPUT_DIR="{{ folder }}"
condor_submit "$jobcard"
{% endfor %}
