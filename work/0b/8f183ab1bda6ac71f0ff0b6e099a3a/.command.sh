#!/bin/bash -ue
fastp \
    -i 2M_PS1269_R1.fastq.gz \
    -o 2M_PS1269_R1.trimmed.fastq.gz \
    -j 2M_PS1269_R1.fastp.json \
    -h 2M_PS1269_R1.fastp.html \
    -w 4
