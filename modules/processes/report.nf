process REPORT {
    tag "$read_id"
    publishDir "${params.outdir}/reports", mode: 'copy'

    cpus 1
    memory '2 GB'

    input:
    tuple val(read_id), path(vcf)
    path table_references

    output:
    path "${read_id}_report.csv", emit: report

    script:
    """
    Rscript ${baseDir}/bin/generate_report.R \\
        --vcf_file ${vcf} \\
        --ref_table ${table_references} \\
        --output_filename ${read_id}_report.csv \\
        --min_read_depth ${params.min_read_depth} \\
        --min_allele_freq ${params.min_allele_freq}
    """
}