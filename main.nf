#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// Results directory
resultsdir = "results/"

// Main workflow
workflow {
    // Input files
    blast_xml = file("data/Trinity_SoF3I50bpF5_paired_mod.blast.xml")

    // Run workflow
    parse_xml(blast_xml)
}

process parse_xml {
    tag "${blast_xml}"
    publishDir "${resultsdir}/tx2gene/",
        mode: "copy"

    input:
    path(blast_xml)

    output:
    path("tx2gene.tsv")

    script:
    """
    parse-xml.py "${blast_xml}" "tx2gene.tsv"
    """
}
