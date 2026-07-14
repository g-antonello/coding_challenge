process ALIGN_MINIMAP2 {
    tag "$read_id"
    publishDir "${params.outdir}/alignments", mode: 'copy'

    cpus 2
    memory '4 GB'

    input:
    tuple val(read_id), path(reads)
    path reference

    output:
    tuple val(read_id), path("${read_id}.sorted.bam"), path("${read_id}.sorted.bam.bai"), emit: bam

    script:
    """
    minimap2 -ax sr -t ${task.cpus} ${reference} ${reads} \\
        | samtools sort -@ ${task.cpus} -o ${read_id}.sorted.bam -

    samtools index ${read_id}.sorted.bam
    """
}

process SAMTOOLS_MPILEUP {
    tag "$read_id"
    publishDir "${params.outdir}/mpileup", mode: 'copy'

    cpus 2
    memory '4 GB'
    
    input:
    tuple val(read_id), path(bam)
    path reference

    output:
    tuple val(read_id), path("${read_id}.mpileup.tsv"), emit: mpileup

    script:
    """
    samtools mpileup \\
        -a \\
        -f ${reference} \\
        # if only the target region is required, leave the rest out
        # -r ${params.target_contig}:${params.target_region} \\
        # filter based on pre-set quality thresholds
        -q ${params.min_map_quality} \\
        -Q ${params.min_base_quality} \\
        ${bam} > ${read_id}.mpileup.tsv
    """

}