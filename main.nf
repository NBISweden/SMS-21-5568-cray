#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// Workflow paths
resultsdir = "results/"

// Main workflow
workflow {
    // References and annotations
    transcript_fasta = file("data/Trinity_SoF3I50bpF5_paired_mod.fasta")
    blast_xml = file("data/Trinity_SoF3I50bpF5_paired_mod.blast.xml")
    mito_fasta = file("data/NC_033509.fasta")

    // Create channels for each read pair across lanes
    Channel
        .fromFilePairs("data/*F2_1/*_L00{1,2}_R1_*.fastq.gz")
        .set{raw_reads_R1}
    Channel
        .fromFilePairs("data/*F2_1/*_L00{1,2}_R2_*.fastq.gz")
        .set{raw_reads_R2}

    // Run the workflow
    parse_xml(blast_xml)
    add_mitochondrion(transcript_fasta,
                      mito_fasta,
                      parse_xml.out.tx2gene_pre_mito)
    salmon_index(add_mitochondrion.out.fasta)
    concatenate_lanes(raw_reads_R1, raw_reads_R2)
}

// Parse the BLAST XML file and get initial tx2gene file
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

// Add mitochondrion genome to FASTA and tx2gene files
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

// Index the transcriptome for use in downstream analyses
process salmon_index {
    publishDir "${resultsdir}/idx/",
        mode: "copy"

    input:
    path(fasta)

    output:
    path("salmon_index", emit: salmon_index)

    script:
    """
    salmon index \
        --transcripts ${fasta} \
        --kmerLen 31 \
        --index "salmon_index"
    """
}

// Concatenate per-sample reads across sequencing lanes
process concatenate_lanes {
    tag "${sample_name}"
    publishDir "${resultsdir}/fastq/",
        mode: "copy"

    input:
    tuple val(sample_name_R1), path(sample_and_lane_fastq_R1)
    tuple val(sample_name_R2), path(sample_and_lane_fastq_R2)

    output:
    tuple val(sample_name), path("*.fastq.gz"), emit: concatenated_reads

    script:
    sample_name = sample_name_R1
    """
    zcat ${sample_and_lane_fastq_R1[0]} ${sample_and_lane_fastq_R1[1]} \
        | gzip > ${sample_name}_R1.fastq.gz
    zcat ${sample_and_lane_fastq_R2[0]} ${sample_and_lane_fastq_R2[1]} \
        | gzip > ${sample_name}_R2.fastq.gz
    """
}
