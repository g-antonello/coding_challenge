#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// use docker/apptainer, best if included in a profile like "cluster/laptop" sort of thing

// Import modules
include { INDEX_REF_DB } from "./modules/processes/ref_db.nf"
include { QC_FASTP } from "./modules/processes/reads_qc.nf"
include { ALIGN_MINIMAP2 } from "./modules/processes/alignments.nf"
include { MPILEUP; VARIANT_CALLING_GATK } from "./modules/processes/variant_calls.nf"
include { REPORT } from "./modules/processes/report.nf"

// other slurm setup


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

	// index the reference database
	INDEX_REF_DB(reference_ch.first())
    
    // 1. QC of reads to have only high quality bases and reads themselves (Phred > = 30)
    QC_FASTP(reads_ch)

    // 2. Align reads using Minimap2 (with Read Groups automatically added)
    ALIGN_MINIMAP2(QC_FASTP.out.trimmed, INDEX_REF_DB.out.indexed_ref)

    // 3. Call variants directly from the BAM files using GATK
    VARIANT_CALLING_GATK(ALIGN_MINIMAP2.out.bam, INDEX_REF_DB.out.indexed_ref)

    
    // 4. Generate a report of the results
    // read reference table before running the report.
    // this is important for report writing, because it contains hard-coded, manually curated positions
    ref_table_ch = Channel.fromPath(params.table_references, checkIfExists: true)
    REPORT(VARIANT_CALLING_GATK.out.vcf, ref_table_ch.first())
	
}