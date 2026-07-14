#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// use docker/apptainer, best if included in a profile like "cluster/laptop" sort of thing

// Import modules
include { FASTP } from "./modules/processes/reads_qc.nf"
include { ALIGN_MINIMAP2; SAMTOOLS_MPILEUP } from "./modules/processes/alignments.nf"
// slurm setup


workflow {
	// initial setup, first gather the reference database
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true)
	// now gather the reads present in the input directory (check nextflow.config for more info)
	// this creates a tuple of read_id and path to the reads file, which is then passed to the FASTP process for quality control
    reads_ch = Channel
        .fromPath(params.reads, checkIfExists: true)
        .map { file ->
            def read_id = file.name.replaceAll(/\.(fastq|fq|fasta|fa|fna)\.gz$/, '')
            tuple(read_id, file)
        }

    FASTP(reads_ch)

    ALIGN_MINIMAP2(FASTP.out.trimmed, reference_ch.first())

	// mpileup generation
	SAMTOOLS_MPILEUP(ALIGN_MINIMAP2.out.bam, reference_ch.first())

}