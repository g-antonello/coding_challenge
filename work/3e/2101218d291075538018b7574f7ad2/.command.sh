#!/bin/bash -ue
fastp \
    -i 2M_PS1269_R2.fastq.gz \
    -o 2M_PS1269_R2.trimmed.fastq.gz \
    -j 2M_PS1269_R2.fastp.json \
    -h 2M_PS1269_R2.fastp.html \
    -w 4
