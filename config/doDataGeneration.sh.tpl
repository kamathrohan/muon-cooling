#!/bin/bash
set -e

jobcard="{{ jobcard }}"
export infile="{{ infile }}"
export ngenerate="{{ ngenerate }}"

{% for folder, outfile in replicas %}
export OUTPUT_DIR="{{ folder }}"
export outfile="{{ outfile }}"
condor_submit "$jobcard"
{% endfor %}
