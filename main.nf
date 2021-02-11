#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// Results directory
resultsdir = "results/"

// Main workflow
workflow {
    // Input files
    transcript_fasta = file("data/Trinity_SoF3I50bpF5_paired_mod.fasta")
    blast_xml = file("data/Trinity_SoF3I50bpF5_paired_mod.blast.xml")
    mito_fasta = file("data/NC_033509.fasta")

    // Run workflow
    parse_xml(blast_xml)
    add_mitochondrion(transcript_fasta,
                      mito_fasta,
                      parse_xml.out.tx2gene_pre_mito)

}

// Process to parse the BLAST XML file and get initial tx2gene file
process parse_xml {

    input:
    path(blast_xml)

    output:
    path("tx2gene-pre-mito.tsv", emit: tx2gene_pre_mito)

    script:
    """
    parse-xml.py ${blast_xml} tx2gene-pre-mito.tsv
    """
}

// Process to add mitochondrion genome to FASTA and tx2gene files
process add_mitochondrion {
    publishDir "${resultsdir}/idx/",
        mode: "copy"

    input:
    path(transcript_fasta)
    path(mito_fasta)
    path(tx2gene_pre_mito)

    output:
    path("pacifastacus-leniusculus.fasta", emit: fasta)
    path("tx2gene.tsv", emit: tx2gene)

    script:
    """
    cat ${transcript_fasta} > pacifastacus-leniusculus.fasta
    cat ${mito_fasta} >> pacifastacus-leniusculus.fasta
    cat ${tx2gene_pre_mito} <(printf "NC_033509\tMT-NC_033509\n") > tx2gene.tsv
    """
}
