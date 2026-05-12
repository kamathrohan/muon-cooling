#!/bin/bash
set -e

jobcard="{{ jobcard }}"
export ngenerate="{{ ngenerate }}"

{% for folder, outfile, skiplines in replicas %}
export OUTPUT_DIR="{{ folder }}"
export outfile="{{ outfile }}"
export skiplines="{{ skiplines }}"
export infile="{{ outfile }}.gmad"
condor_submit "$jobcard"
{% endfor %}
