#!/usr/bin/env Rscript

# TMM-normalized donor pseudobulks after excluding the two subtypes that have
# fewer than five recovered cells in at least one BPD donor. This sensitivity
# addresses whether the all-endothelial score is driven by abCap or systemic
# venous EC composition.

script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
script_path <- normalizePath(
  sub("^--file=", "", script_argument[[1]]),
  winslash = "/",
  mustWork = TRUE
)
project_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  winslash = "/",
  mustWork = TRUE
)
.libPaths(c(file.path(project_root, ".R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

processed_dir <- file.path(project_root, "04_processed", "GSE275938")
result_dir <- file.path(
  project_root,
  "07_results",
  "human_endothelial_subtypes"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

counts <- fread(
  file.path(
    processed_dir,
    "GSE275938_endothelial_subtype_pseudobulk_counts.tsv.gz"
  )
)
metadata <- fread(
  file.path(
    processed_dir,
    "GSE275938_endothelial_subtype_pseudobulk_metadata.tsv"
  )
)

sample_order <- c("BPD 1", "BPD 2", "BPD+PH 1", "BPD+PH 2")
common_subtypes <- c(
  "gCap",
  "aCap",
  "Arterial EC",
  "Pulmonary venous EC",
  "Lymphatic"
)

counts_matrix <- as.matrix(counts[, -"pseudobulk_id"])
rownames(counts_matrix) <- counts$pseudobulk_id
storage.mode(counts_matrix) <- "integer"
genes <- colnames(counts_matrix)

aggregated <- matrix(
  0L,
  nrow = length(genes),
  ncol = length(sample_order),
  dimnames = list(genes, sample_order)
)
metadata_rows <- list()

for (sample_name in sample_order) {
  selected <- metadata[
    sample == sample_name & subtype %in% common_subtypes
  ]
  if (nrow(selected) != length(common_subtypes)) {
    stop("Incomplete common-subtype set for ", sample_name)
  }
  donor_counts <- colSums(
    counts_matrix[selected$pseudobulk_id, , drop = FALSE]
  )
  aggregated[, sample_name] <- donor_counts
  metadata_rows[[sample_name]] <- data.table(
    sample = sample_name,
    condition = selected$condition[[1]],
    included_subtypes = paste(common_subtypes, collapse = ";"),
    excluded_subtypes = "abCap;Systemic venous EC",
    cells = sum(selected$cells),
    library_size = sum(donor_counts),
    genes_detected = sum(donor_counts > 0)
  )
}

y <- DGEList(
  counts = aggregated,
  group = factor(
    c("BPD", "BPD", "BPD+PH", "BPD+PH"),
    levels = c("BPD", "BPD+PH")
  )
)
y <- calcNormFactors(y, method = "TMM")
log2cpm <- cpm(
  y,
  log = TRUE,
  prior.count = 0.5,
  normalized.lib.sizes = TRUE
)

metadata_output <- rbindlist(metadata_rows)
metadata_output[, norm_factor := y$samples$norm.factors]
metadata_output[, effective_library_size := y$samples$lib.size * norm_factor]

raw_wide <- as.data.table(t(aggregated), keep.rownames = "sample")
log_wide <- as.data.table(t(log2cpm), keep.rownames = "sample")
raw_long <- melt(
  raw_wide,
  id.vars = "sample",
  variable.name = "gene",
  value.name = "raw_count"
)
log_long <- melt(
  log_wide,
  id.vars = "sample",
  variable.name = "gene",
  value.name = "log2CPM_TMM"
)
expression_output <- merge(
  raw_long,
  log_long,
  by = c("sample", "gene"),
  sort = FALSE
)
expression_output <- merge(
  expression_output,
  metadata_output[, .(sample, condition)],
  by = "sample",
  sort = FALSE
)
expression_output[, aggregate := "Common five endothelial subtypes"]
expression_output[, detected := raw_count > 0]
setcolorder(
  expression_output,
  c(
    "sample",
    "condition",
    "aggregate",
    "gene",
    "raw_count",
    "detected",
    "log2CPM_TMM"
  )
)

fwrite(
  metadata_output,
  file.path(
    result_dir,
    "GSE275938_common_five_subtypes_aggregate_metadata.tsv"
  ),
  sep = "\t",
  na = "NA"
)
fwrite(
  expression_output,
  file.path(
    result_dir,
    "GSE275938_common_five_subtypes_aggregate_expression.tsv.gz"
  ),
  sep = "\t",
  na = "NA"
)
writeLines(
  capture.output(sessionInfo()),
  file.path(
    result_dir,
    "GSE275938_common_five_subtypes_TMM_sessionInfo.txt"
  )
)

print(metadata_output)
