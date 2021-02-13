#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// Workflow paths
resultsdir = "results/"

// Main workflow
workflow {
    // References and annotations
    transcript_fasta = file("data/Trinity_SoF3I50bpF5_paired_mod.fasta")
    blast_xml = file("data/Trinity_SoF3I50bpF5_paired_mod.blast.xml")
    mitochondrion_fasta = file("data/NC_033509.fasta")

    // Create channels for each read pair across lanes
    Channel
        .fromFilePairs("data/*/*_L00{1,2}_R1_*.fastq.gz")
        .set{raw_reads_R1}
    Channel
        .fromFilePairs("data/*/*_L00{1,2}_R2_*.fastq.gz")
        .set{raw_reads_R2}

    // Run the workflow
    build_tx2gene_and_fasta(blast_xml,
                  transcript_fasta,
                  mitochondrion_fasta)
    salmon_index(build_tx2gene_and_fasta.out.fasta)
    concatenate_lanes(raw_reads_R1,
                      raw_reads_R2)
    quantify_expression(concatenate_lanes.out.concatenated_reads,
                        salmon_index.out.salmon_index,
                        build_tx2gene_and_fasta.out.tx2gene)
}

// Build the tx2gene mapping and transcriptome + mitochondrion FASTA
process build_tx2gene_and_fasta {
    input:
    path(blast_xml)
    path(transcript_fasta)
    path(mitochondrion_fasta)

    output:
    path("tx2gene.tsv", emit: tx2gene)
    path("pacifastacus-leniusculus.fasta", emit: fasta)

    script:
    """
    # Parse the BLAST XML file and grab existing annotations
    parse-xml.py ${blast_xml} tx2gene-annotated.tsv

    # Add transcripts without gene annotations
    join -a1 -a2 -e- \
        <(sort -k1,1 tx2gene-annotated.tsv) \
        <(grep ">" ${transcript_fasta} | sed 's/>//g' | sort -k1,1) \
        > tx2gene-joined.tsv

    # Add transcript IDs as gene IDs when no annotation is available
    paste \
        <(cut -d ' ' -f 1 tx2gene-joined.tsv) \
        <(cut -d ' ' -f 2 tx2gene-joined.tsv) \
        > tx2gene-pre-mito.tsv

    # Add mitochondrion genome entry to tx2gene
    cat tx2gene-pre-mito.tsv <(printf "NC_033509.1\tMT-NC_033509\n") \
        > tx2gene.tsv

    # Concatenate transcriptome and mitochondrion genome
    cat ${transcript_fasta} > pacifastacus-leniusculus.fasta
    cat ${mitochondrion_fasta} >> pacifastacus-leniusculus.fasta
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

// Quantify the raw reads with Alevin
process quantify_expression {
    tag "${sample_name}"
    publishDir "${resultsdir}/expression/",
        mode: "copy"

    input:
    tuple val(sample_name), path(sample_fastq)
    path(salmon_index)
    path(tx2gene)

    output:
    path("${sample_name}", emit: alevin_output)

    script:
    """
    salmon alevin \
        --libType ISR \
        --mates1 ${sample_fastq[0]} \
        --mates2 ${sample_fastq[1]} \
        --chromiumV3 \
        --index ${salmon_index} \
        --output ${sample_name} \
        --tgMap ${tx2gene}
    """
}
