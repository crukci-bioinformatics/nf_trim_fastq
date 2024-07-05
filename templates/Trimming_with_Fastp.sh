#!/bin/bash

# This script is used to trim adapters from fastq files using fastp

set -exu 

echo Start `date`

fastp \
    -i ${fastqR1In} \
    -I ${fastqR2In} \
    -o ${fastqR1Out} \
    -O ${fastqR2Out} \
    --detect_adapter_for_pe \
    -g \
    -x \
    -j ${fastpJson} \
    -h ${fastpHtml} \
    --thread 8 \
    ${params.fastpOptions}

# some fastp options:
# -A - Disable adapter trimming. Adapter trimming is enabled by default.
# -x - enable polyX trimming in 3' ends.
# -g - enable polyG tail trimming
# -G - disable polyG tail trimming, by default trimming is automatically enabled
#      for Illumina NextSeq/NovaSeq data (if detected)
# -Q - disable_quality_filtering, quality filtering is enabled by default.
# -L - disable length filtering, by default filtering is enabled (with a
#      minminal length of 15bp)
# -y - enable low complexity filtering, by default filtering is disabled

fastp --version