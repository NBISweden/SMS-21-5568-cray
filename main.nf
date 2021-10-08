#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// Workflow paths
resultsdir = "results/"

workflow {
    // Analysis workflow

    // NCBI transcriptome
    ncbi_transcriptome = file("https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GB/YW/GBYW01/GBYW01.1.fsa_nt.gz")

    // Other data files
    cell_type_markers = file("data/cell-type-markers.csv")
    features_to_plot = file("data/features-to-plot.tsv")
    cluster_cell_types = file("data/cluster-cell-types.tsv")

    // Report files
    report_1_qc = file("bin/report-1-quality-controls.Rmd")
    report_2_integration = file("bin/report-2-integration-and-markers.Rmd")
    report_3_de = file("bin/report-3-differential-expression.Rmd")
    report_4_plots = file("bin/report-4-features-and-celltypes.Rmd")
    report_5_ti = file("bin/report-5-trajectory-inference.Rmd")

    // Input files
    Channel
        .fromPath("results/expression/TJ-2700-*/", type: "dir")
        .toList()
        .set{gene_expression}

    // Run the workflow
    map_transcriptome_ids(ncbi_transcriptome)
    merge_cells(gene_expression)
    quality_controls(report_1_qc, merge_cells.out.seurat_object)
    integration(report_2_integration,
                quality_controls.out.seurat_object,
                map_transcriptome_ids.out.transcriptome_id_map,
                cell_type_markers)
    differential_expression(report_3_de,
                            integration.out.seurat_object)
    features_and_celltypes(report_4_plots,
                           integration.out.seurat_object,
                           features_to_plot,
                           cluster_cell_types)
    trajectory_inference(report_5_ti,
                         features_and_celltypes.out.seurat_object)
}

process map_transcriptome_ids {
    // Map NCBI <-> internal transcriptome IDs
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
    gzip --decompress < ${ncbi_transcriptome} > tmp.fasta
    grep ">" tmp.fasta \
        | cut -d ' ' -f 1,5 \
        | tr -d '>' \
        | tr ' ' '\t' \
        | cut -d '_' -f 1,2 \
        >> transcriptome_id_map.tsv
    """
}

process merge_cells {
    // Merge same-sample cells
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

process quality_controls {
    // Basic QC plots and filtering
    publishDir "${resultsdir}",
        mode: "copy",
        saveAs: { filename ->
            filename.indexOf(".html") > 0 ? \
                "${filename}" : "seurat/01-quality_controls/${filename}"
        }

    input:
    path(report)
    path(seurat_object)

    output:
    path("report-1-quality-controls.html")
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

process integration {
    // Sample integration and known cell type marker expression
    publishDir "${resultsdir}",
        mode: "copy",
        saveAs: { filename ->
            filename.indexOf(".html") > 0 ? \
                "${filename}" : "seurat/02-integration-and-markers/${filename}"
        }

    input:
    path(report)
    path(seurat_object)
    path(transcriptome_id_map)
    path(cell_type_markers)

    output:
    path("report-2-integration-and-markers.html")
    path("seurat-integration.rds", emit: seurat_object)

    script:
    """
    #!/usr/bin/env Rscript
    parameters <- list(root_directory       = getwd(),
                       input_seurat_object  = "${seurat_object}",
                       output_seurat_object = "seurat-integration.rds",
                       transcriptome_id_map = "${transcriptome_id_map}",
                       cell_type_markers    = "${cell_type_markers}")
    output_report <- gsub(".Rmd", ".html", basename("${report}"))
    rmarkdown::render("${report}",
                      params      = parameters,
                      output_dir  = getwd(),
                      output_file = output_report)
    """
}

process differential_expression {
    // Differential expression analyses
    publishDir "${resultsdir}",
        mode: "copy",
        saveAs: { filename ->
            filename.indexOf(".html") > 0 ? \
                "${filename}" : "seurat/03-differential-expression/${filename}"
        }

    input:
    path(report)
    path(seurat_object)

    output:
    path("report-3-differential-expression.html")
    path("*.tsv", optional: true)

    script:
    """
    #!/usr/bin/env Rscript
    parameters <- list(root_directory       = getwd(),
                       input_seurat_object  = "${seurat_object}")
    output_report <- gsub(".Rmd", ".html", basename("${report}"))
    rmarkdown::render("${report}",
                      params      = parameters,
                      output_dir  = getwd(),
                      output_file = output_report)
    """
}

process features_and_celltypes {
    publishDir "${resultsdir}/",
        mode: "copy",
        saveAs: { filename ->
            filename.indexOf(".html") > 0 ? \
                "${filename}" : "seurat/04-features-and-celltypes/${filename}"
        }

    input:
    path(report)
    path(seurat_object)
    path(features_to_plot)
    path(cluster_cell_types)

    output:
    path("report-4-features-and-celltypes.html")
    path("seurat-features-and-celltypes.rds", emit: seurat_object)

    script:
    """
    #!/usr/bin/env Rscript
    parameters <- list(root_directory       = getwd(),
                       input_seurat_object  = "${seurat_object}",
                       output_seurat_object = "seurat-features-and-celltypes.rds",
                       features_to_plot     = "${features_to_plot}",
                       cluster_cell_types   = "${cluster_cell_types}")
    output_report <- gsub(".Rmd", ".html", basename("${report}"))
    rmarkdown::render("${report}",
                      params      = parameters,
                      output_dir  = getwd(),
                      output_file = output_report)
    """
}

process trajectory_inference {
    container "nbis-5568-tradeseq"
    publishDir "${resultsdir}/",
        mode: "copy",
        saveAs: { filename ->
            filename.indexOf(".html") > 0 ? \
                "${filename}" : "seurat/05-trajectory-inference/${filename}"
        }

    input:
    path(report)
    path(seurat_object)

    output:
    path("report-5-trajectory-inference.html")
    // path("seurat-trajectory-inference.rds", emit: seurat_object)

    script:
    """
    #!/usr/bin/env Rscript
    parameters <- list(root_directory       = getwd(),
                       input_seurat_object  = "${seurat_object}")
    output_report <- gsub(".Rmd", ".html", basename("${report}"))
    rmarkdown::render("${report}",
                      params      = parameters,
                      output_dir  = getwd(),
                      output_file = output_report)
    """
}

workflow pre_processing {
    // Pre-processing workflow

    // References and annotations
    blast_xml = file("data/Trinity_SoF3I50bpF5_paired_mod.blast.xml")
    transcript_fasta = file("data/Trinity_SoF3I50bpF5_paired_mod.fasta")
    manual_tx2gene = file("data/manual-tx2gene.tsv")
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
                            manual_tx2gene,
                            mitochondrion_fasta)
    salmon_index(build_tx2gene_and_fasta.out.fasta)
    concatenate_lanes(raw_reads_R1,
                      raw_reads_R2)
    quantify_expression(concatenate_lanes.out.concatenated_reads,
                        salmon_index.out.salmon_index,
                        build_tx2gene_and_fasta.out.tx2gene)
    alevin_quality_controls(quantify_expression.out.alevin_output)
}

process build_tx2gene_and_fasta {
    // Build the tx2gene mapping and transcriptome + mitochondrion FASTA
    publishDir "${resultsdir}/idx/",
        mode: "copy"

    input:
    path(blast_xml)
    path(transcript_fasta)
    path(manual_tx2gene)
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
    # (remove `_seqZ` part for summing to "gene" level)
    paste \
        <(cut -d ' ' -f 1 tx2gene-joined.tsv) \
        <(cut -d ' ' -f 2 tx2gene-joined.tsv | cut -f 1,2 -d '_') \
        > tx2gene-pre-manual.tsv

    # Add manual mappings from the group
    join -a1 -a2 -e- <(sort ${manual_tx2gene}) <(sort tx2gene-pre-manual.tsv) \
        | tr ' ' '\t' \
        | cut -f 1,2 \
        > tx2gene-pre-mito.tsv

    # Add mitochondrion genome entry to tx2gene
    cat tx2gene-pre-mito.tsv <(printf "NC_033509.1\tMT-NC_033509\n") \
        > tx2gene.tsv

    # Concatenate transcriptome and mitochondrion genome
    cat ${transcript_fasta} > pacifastacus-leniusculus.fasta
    cat ${mitochondrion_fasta} >> pacifastacus-leniusculus.fasta
    """
}

process salmon_index {
    // Index the transcriptome for use in downstream analyses
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

process concatenate_lanes {
    // Concatenate per-sample reads across sequencing lanes
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

process quantify_expression {
    // Quantify the raw reads with Alevin
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

process alevin_quality_controls {
    // Perform basic quality controls of the Alevin gene expression data
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
