process FASTP {
    tag "$read_id"
    publishDir "${params.outdir}/qc/fastp", mode: 'copy'

    cpus 4
    memory '8 GB'

    input:
    tuple val(read_id), path(reads)

    output:
    tuple val(read_id), path("${read_id}.trimmed.fastq.gz"), emit: trimmed
    path "${read_id}.fastp.json", emit: json
    path "${read_id}.fastp.html", emit: html

    script:
    """
    fastp \\
        -i ${reads} \\
        -o ${read_id}.trimmed.fastq.gz \\
        -j ${read_id}.fastp.json \\
        -h ${read_id}.fastp.html \\
        -w ${task.cpus}
    """
}
