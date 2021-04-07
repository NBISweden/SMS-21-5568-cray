#!/usr/bin/env Rscript

# Argument parser
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("input_files",
                    type    = "character",
                    help    = "Gene expression input files")
parser$add_argument("output_path",
                    type    = "character",
                    help    = "Output path for the final Seurat object")
args <- parser$parse_args()

# Load packages
suppressPackageStartupMessages({
    library("tximport")
    library("Seurat")
})

# Function for creating Seurat objects
create_seurat_object <- function(file) {

    # Get sample name from the file name
    sample_name <- strsplit(file, "/")[[1]][3]

    # Import data
    txi <- tximport(file, type = "alevin")

    # Create the Seurat object
    s_obj <- CreateSeuratObject(counts       = txi$counts,
                                min.cells    = 0,
                                min.features = 0,
                                project      = sample_name)
    return(s_obj)
}

# Function for merging same-sample cells
merge_cells <- function(objects, sample_ids, final_sample_id) {

    # Merge objects
    last_three_objects <- c(objects[[2]], objects[[3]], objects[[4]])
    s_obj <- merge(objects[[1]],
                   y            = last_three_objects,
                   add.cell.ids = sample_ids)

    # Get cell IDs
    cell_ids <- sub(".*_", "", colnames(s_obj))

    # Sum expression across same cells
    sums <- rowsum(as.matrix(Matrix::t(s_obj@assays$RNA@counts)), cell_ids)
    sums <- Matrix::Matrix(sums, sparse = TRUE)
    sums <- Matrix::t(sums)

    # Create new Seurat object from summed expression
    s_obj_all <- CreateSeuratObject(counts       = sums,
                                    min.cells    = 3,
                                    min.features = 10,
                                    project      = final_sample_id)
}

# List all expression files
files <- sort(strsplit(args$input_files, " ")[[1]])
files <- file.path(files, "alevin/quants_mat.gz")

# Create individual Seurat objects (s1-4)
s1 <- create_seurat_object(files[1])
s2 <- create_seurat_object(files[2])
s3 <- create_seurat_object(files[3])
s4 <- create_seurat_object(files[4])

# Merge same-sample cells
s_obj_sample_1 <- merge_cells(list(s1, s2, s3, s4),
                              c("S1", "S2", "S3", "S4"),
                              "S1")
rm(s1, s2, s3, s4)
invisible(gc())

# Create individual Seurat objects (s5-8)
s5 <- create_seurat_object(files[5])
s6 <- create_seurat_object(files[6])
s7 <- create_seurat_object(files[7])
s8 <- create_seurat_object(files[8])
s_obj_sample_2 <- merge_cells(list(s5, s6, s7, s8),
                              c("S5", "S6", "S7", "S8"),
                              "S2")
rm(s5, s6, s7, s8)
invisible(gc())

# Merge per-sample Seurat objects
s_obj <- merge(x            = s_obj_sample_1,
               y            = s_obj_sample_2,
               add.cell.ids = c("S1", "S2"))
rm(s_obj_sample_1, s_obj_sample_2)
invisible(gc())

# Save the final Seurat object
saveRDS(s_obj, file = args$output_path)
