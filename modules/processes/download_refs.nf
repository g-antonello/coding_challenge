// read accessions from file

process download_refs (accession){
    input: ${accessions}

    command:
    """

    datasets download genome accession ${accession} --include genome,gff3

    // unzip 
    unzip ncbi_dataset.zip

    // fing 23S ribosomal RNA sequences in gff3 
    grep -i "23S ribosomal RNA" ncbi_dataset/data/${accession}/genomic.gff | grep "type="rRNA" > ${accession}_23S_header.txt

    seqkit subseq --gtf genomic.gff --feature-type rRNA genomic.fna > all_rrna.fa

    // search "23S ribosomal RNA"
    grep "23S ribosomal RNA" ncbi_dataset/data/${accession}/genomic.gff

    // clean up

    """
}

refs_list = file('')

// final output should be a path to a reference that can be used for alignment.
params.reference = "data/ref_23S_database.fna"