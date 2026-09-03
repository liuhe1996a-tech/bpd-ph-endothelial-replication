#!/usr/bin/env Rscript

# Export the independent GSE243129 WT capillary cells into a language-neutral
# sparse matrix for the R9 perturbation-prediction benchmark.  The deposited
# object is gzip-compressed twice; only the project-local temporary copy is
# decompressed and it is removed when the script exits.

suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
})

options(stringsAsFactors = FALSE, width = 140)

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
root <- if (length(args)) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else
  normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

raw_file <- file.path(root, "raw", "GSE243129", "GSE243129_combined_labeled.RDS.gz")
output_dir <- file.path(root, "raw", "GSE243129_virtual_cell")
tmp_dir <- file.path(root, "tmp")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(raw_file))

copy_connection <- function(input, output, chunk = 16L * 1024L * 1024L) {
  src <- gzfile(input, "rb")
  dst <- file(output, "wb")
  on.exit({
    try(close(src), silent = TRUE)
    try(close(dst), silent = TRUE)
  }, add = TRUE)
  repeat {
    bytes <- readBin(src, what = "raw", n = chunk)
    if (!length(bytes)) break
    writeBin(bytes, dst)
  }
  invisible(output)
}

gzip_file <- function(input, output, chunk = 16L * 1024L * 1024L) {
  src <- file(input, "rb")
  dst <- gzfile(output, "wb", compression = 6)
  on.exit({
    try(close(src), silent = TRUE)
    try(close(dst), silent = TRUE)
  }, add = TRUE)
  repeat {
    bytes <- readBin(src, what = "raw", n = chunk)
    if (!length(bytes)) break
    writeBin(bytes, dst)
  }
  invisible(output)
}

inner_file <- file.path(tmp_dir, "GSE243129_virtual_cell_inner.RDS.gz")
matrix_file <- file.path(tmp_dir, "GSE243129_WT_capillary_counts_cells_by_genes.mtx")
copy_connection(raw_file, inner_file)
on.exit(unlink(c(inner_file, matrix_file), force = TRUE), add = TRUE)

object <- readRDS(inner_file)
stopifnot(inherits(object, "Seurat"), "RNA" %in% names(object@assays))
meta <- object@meta.data
meta$cell_id <- rownames(meta)
meta$deposited_cell_type <- as.character(object@active.ident)
required <- c("cell_id", "deposited_cell_type", "Ex", "Time", "Type", "Rep", "Lane")
if (!all(required %in% colnames(meta))) {
  stop("Missing metadata fields: ", paste(setdiff(required, colnames(meta)), collapse = ", "))
}

meta$Age <- as.character(meta$Time)
meta$Oxygen <- ifelse(meta$Ex == "O2", "Hyperoxia", "Normoxia")
meta$genotype <- ifelse(meta$Type == "CTRL", "WT", "Tgfbr2_KO")
meta$animal_id <- paste(meta$Time, meta$Type, meta$Ex, meta$Rep, sep = "_")
meta$CellType <- ifelse(
  meta$deposited_cell_type == "gCap", "Cap",
  ifelse(meta$deposited_cell_type == "aCap", "Cap-a", NA_character_)
)

keep <- meta$genotype == "WT" & !is.na(meta$CellType) & meta$Age %in% c("P7", "P14")
if (!any(keep)) stop("No WT capillary cells were found")
export_meta <- meta[keep, c(
  "cell_id", "animal_id", "Age", "Oxygen", "CellType", "deposited_cell_type",
  "genotype", "Rep", "Lane"
), drop = FALSE]
counts <- GetAssayData(object, assay = "RNA", layer = "counts")
counts <- counts[, export_meta$cell_id, drop = FALSE]
stopifnot(identical(colnames(counts), export_meta$cell_id))
if (any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-8)) {
  stop("The deposited count layer is not non-negative integer-valued")
}

writeMM(t(counts), matrix_file)
gzip_file(
  matrix_file,
  file.path(output_dir, "GSE243129_WT_capillary_counts_cells_by_genes.mtx.gz")
)
write.table(
  data.frame(gene = rownames(counts)),
  file.path(output_dir, "GSE243129_WT_capillary_genes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
metadata_connection <- gzfile(
  file.path(output_dir, "GSE243129_WT_capillary_cell_metadata.tsv.gz"), "wt"
)
write.table(
  export_meta,
  metadata_connection,
  sep = "\t", quote = FALSE, row.names = FALSE
)
close(metadata_connection)

animal_audit <- aggregate(
  cell_id ~ animal_id + Age + Oxygen + CellType,
  data = export_meta,
  FUN = length
)
names(animal_audit)[names(animal_audit) == "cell_id"] <- "n_cells"
animal_audit <- animal_audit[order(
  animal_audit$Age, animal_audit$Oxygen, animal_audit$CellType, animal_audit$animal_id
), ]
write.table(
  animal_audit,
  file.path(output_dir, "GSE243129_WT_capillary_animal_cell_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

audit <- data.frame(
  item = c(
    "accession", "model_role", "n_cells", "n_genes", "n_animals", "ages",
    "oxygen_groups", "cell_types", "biological_unit", "lane_handling"
  ),
  value = c(
    "GSE243129", "independent external perturbation-prediction cohort",
    ncol(counts), nrow(counts), length(unique(export_meta$animal_id)),
    paste(sort(unique(export_meta$Age)), collapse = ";"),
    paste(sort(unique(export_meta$Oxygen)), collapse = ";"),
    paste(sort(unique(export_meta$CellType)), collapse = ";"),
    "mouse", "two sequencing lanes merged through the deposited Rep identifier"
  )
)
write.table(
  audit,
  file.path(output_dir, "GSE243129_WT_capillary_export_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Exported", ncol(counts), "cells,", nrow(counts), "genes and",
    length(unique(export_meta$animal_id)), "animals\n")
print(animal_audit)
