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