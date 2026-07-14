#!/bin/bash -ue
minimap2 -ax sr -t 2 ref_23S_database.fna 2M_PS1269_R1.trimmed.fastq.gz \
    | samtools sort -@ 2 -o 2M_PS1269_R1.sorted.bam -

samtools index 2M_PS1269_R1.sorted.bam
