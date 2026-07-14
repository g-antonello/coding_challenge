process MPILEUP {
    tag "$read_id"
    publishDir "${params.outdir}/mpileup", mode: 'copy'

    cpus 2
    memory '4 GB'
    
    input:
    tuple val(read_id), path(bam), path(bai)
    path reference

    output:
    tuple val(read_id), path("${read_id}.mpileup.tsv"), emit: mpileup

    /* output .tsv files are 6-column .tsv
    * Fields are: chrom, pos, ref_base, depth, read_bases, base_qualities
    * The read-bases column (./, = matches ref forward/reverse; 
    * ACGTacgt = mismatch
    * ^]/$ = read start/end markers; * = deletion)
    */

    script:

    
    """
    samtools mpileup \
        -f ${reference} \
        -q ${params.min_map_quality} \
        -Q ${params.min_base_quality} \
        ${bam} > ${read_id}.mpileup.tsv
    """

}

process VARIANT_CALLING_GATK {
    tag "$read_id"
    publishDir "${params.outdir}/variants", mode: 'copy'

    // GATK can be resource-intensive; increased memory to 4 GB is safer
    cpus 2
    memory '4 GB'

    input:
    // 1. Accept the BAM and BAI (index) from ALIGN_MINIMAP2
    tuple val(read_id), path(bam), path(bai)
    // 2. Accept the 3-file indexed reference bundle
    tuple path(ref_fna), path(ref_fai), path(ref_dict)

    output:
    // Output a standard GATK VCF file
    tuple val(read_id), path("${read_id}.vcf"), emit: vcf

    // mutect2 is considered the good way to go for bacterial variant calling
    script:
    """
    gatk Mutect2 \\
        -R ${ref_fna} \\
        -I ${bam} \\
        -O ${read_id}.vcf
    """
}

// below older implementations
process VARIANT_CALLING {
    tag "$read_id"
    publishDir "${params.outdir}/variants", mode: 'copy'

    cpus 1
    memory '2 GB'

    input:
    tuple val(read_id), path(mpileup_tsv)

    output:
    path "${read_id}.variant_calls.tsv", emit: calls

    script:
    """
    python3 ${projectDir}/bin/call_resistance_mutations.py \\
        --pileup ${mpileup_tsv} \\
        --sample-id ${read_id} \\
        --min-coverage ${params.min_coverage} \\
        --min-variant-freq ${params.min_variant_freq} \\
        --output ${read_id}.variant_calls.tsv
    """
}

// this below is a process that merges all the variant calls from different samples into a single file, which can be useful for downstream analysis or reporting.
// not implemented yet, but can be added later if needed.

process MERGE_VARIANT_CALLS {
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    path(call_files)

    output:
    path "all_samples.variant_calls.tsv"

    script:
    """
    head -n 1 ${call_files[0]} > all_samples.variant_calls.tsv
    for f in ${call_files}; do
        tail -n +2 \$f >> all_samples.variant_calls.tsv
    done
    """
}