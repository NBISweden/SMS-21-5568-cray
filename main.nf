#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// Workflow paths
resultsdir = "results/"

// Analysis workflow
workflow {
    // Report files
    report_1_qc = file("bin/report_1_quality_controls.Rmd")
    report_2_de = file("bin/report_2_differential_expression.Rmd")
    ncbi_transcriptome = file("https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GB/YW/GBYW01/GBYW01.1.fsa_nt.gz")

    // Input files
    Channel
        .fromPath("results/expression/TJ-2700-*/", type: "dir")
        .toList()
        .set{gene_expression}

    // Run the workflow
    map_transcriptome_ids(ncbi_transcriptome)
    merge_cells(gene_expression)
    quality_controls(report_1_qc, merge_cells.out.seurat_object)
    differential_expression(report_2_de, quality_controls.out.seurat_object)
}

// Map NCBI <-> internal transcriptome IDs
process map_transcriptome_ids {
    publishDir "${resultsdir}/idx/",
        mode: "copy"

    input:
    path(ncbi_transcriptome)

    output:
    path("transcriptome_id_map.tsv", emit: transcriptome_id_map)

    script:
    """
    # Initialise mapping file
    printf "ncbi_id\toriginal_id\n" \
        > transcriptome_id_map.tsv

    # Get ID mappings
    gzcat GBYW01.1.fsa_nt.gz > tmp.fasta
    grep ">" tmp.fasta \
        | cut -d ' ' -f 1,5 \
        | tr -d '>' \
        | tr ' ' '\t' \
        >> transcriptome_id_map.tsv
    """
}

// Merge same-sample cells
process merge_cells {
    publishDir "${resultsdir}/seurat/00-merged-cells",
        mode: "copy"

    input:
    path(gene_expression)

    output:
    path("seurat-merged.rds", emit: seurat_object)

    script:
    """
    merge_cells.R "${gene_expression}" "seurat-merged.rds"
    """
}

// Basic QC plots and filtering
process quality_controls {
    publishDir "${resultsdir}/seurat/01-quality_controls/",
        mode: "copy"

    input:
    path(report)
    path(seurat_object)

    output:
    path("report_1_quality_controls.html")
    path("seurat-qc.rds", emit: seurat_object)

    script:
    """
    #!/usr/bin/env Rscript
    parameters <- list(root_directory       = getwd(),
                       input_seurat_object  = "${seurat_object}",
                       output_seurat_object = "seurat-qc.rds")
    output_report <- gsub(".Rmd", ".html", basename("${report}"))
    rmarkdown::render("${report}",
                      params      = parameters,
                      output_dir  = getwd(),
                      output_file = output_report)
    """
}

process differential_expression {
    publishDir "${resultsdir}/seurat/02-differential-expression",
        mode: "copy"

    input:
    path(report)
    path(seurat_object)

    output:
    path("report_2_differential_expression.html")
    path("seurat-de.rds")

    script:
    """
    #!/usr/bin/env Rscript
    parameters <- list(root_directory       = getwd(),
                       input_seurat_object  = "${seurat_object}",
                       output_seurat_object = "seurat-de.rds")
    output_report <- gsub(".Rmd", ".html", basename("${report}"))
    rmarkdown::render("${report}",
                      params      = parameters,
                      output_dir  = getwd(),
                      output_file = output_report)
    """
}

// Pre-processing workflow
workflow pre_processing {
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
    alevin_quality_controls(quantify_expression.out.alevin_output)
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
    tuple val(sample_name), path("${sample_name}"), emit: alevin_output

    script:
    """
    salmon alevin \
        --libType ISR \
        --mates1 ${sample_fastq[0]} \
        --mates2 ${sample_fastq[1]} \
        --chromiumV3 \
        --index ${salmon_index} \
        --output ${sample_name} \
        --tgMap ${tx2gene} \
        --dumpFeatures
    """
}

// Perform basic quality controls of the Alevin gene expression data
process alevin_quality_controls {
    tag "${sample_name}"
    publishDir "${resultsdir}/expression/${sample_name}",
        mode: "copy"

    input:
    tuple val(sample_name), path(alevin_output)

    output:
    path("alevinReport.${sample_name}.html")

    script:
    """
    #!/usr/bin/env Rscript
    library("alevinQC")
    alevinQCReport(baseDir        = "${sample_name}",
                   outputDir      = "./",
                   ignorePandoc   = TRUE,
                   forceOverwrite = TRUE,
                   outputFormat   = "html_document",
                   sampleId       = "${sample_name}",
                   outputFile     = "alevinReport.${sample_name}.html")
    """
}
