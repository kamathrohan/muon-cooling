#!/bin/bash
set -e

jobcard="{{ jobcard }}"
export ngenerate="{{ ngenerate }}"

{% for folder, outfile in replicas %}
export OUTPUT_DIR="{{ folder }}"
export outfile="{{ outfile }}"
export infile="{{ outfile }}.gmad"
condor_submit "$jobcard"
{% endfor %}
