process ALIGN_MINIMAP2 {
    tag "$read_id"
    publishDir "${params.outdir}/alignments", mode: 'copy'

    cpus 2
    memory '4 GB'

    input:
    tuple val(read_id), path(reads)
    // input from the indexing is a 3-file bundle: the fasta, the fai, and the dict
    // we only need the fasta, but let's be complete and input everything, the select later
    tuple path(ref_fna), path(ref_fai), path(ref_dict)

    output:
    tuple val(read_id), path("${read_id}.sorted.bam"), path("${read_id}.sorted.bam.bai"), emit: bam

    script:
    """
    minimap2 -ax sr -R '@RG\\tID:${read_id}\\tLB:lib1\\tPL:illumina\\tPU:unit1\\tSM:${read_id}' -t ${task.cpus} ${ref_fna} ${reads} \\
        | samtools sort -@ ${task.cpus} -o ${read_id}.sorted.bam -
    
    samtools index ${read_id}.sorted.bam
    """
}