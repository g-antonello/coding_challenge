// read accessions from file
process INDEX_REF_DB {
    input:
    path fasta

    output:
    // This bundles all three files into a single channel
    tuple path(fasta), path("${fasta}.fai"), path("${fasta.baseName}.dict"), emit: indexed_ref

    script:
    """
    samtools faidx ${fasta}
    gatk CreateSequenceDictionary -R ${fasta}
    """
    
}