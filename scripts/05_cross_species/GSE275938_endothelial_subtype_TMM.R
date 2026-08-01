options(stringsAsFactors = FALSE)

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
  Sys.getenv(
    "BPD_PH_PROJECT_ROOT",
    unset = file.path(dirname(script_path), "..", "..")
  ),
  winslash = "/",
  mustWork = TRUE
)
local_library <- file.path(project_root, ".R_library")
.libPaths(c(local_library, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

processed_dir <- file.path(project_root, "04_processed", "GSE275938")
result_dir <- file.path(
  project_root, "07_results", "human_endothelial_subtypes"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

count_table <- fread(
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
locked <- fread(
  file.path(
    project_root,
    "07_results",
    "cross_species",
    "locked_mouse_signature_human_BPD_PH_evidence.tsv"
  )
)

pseudobulk_ids <- count_table$pseudobulk_id
gene_names <- setdiff(names(count_table), "pseudobulk_id")
count_matrix <- t(as.matrix(count_table[, ..gene_names]))
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- gene_names
colnames(count_matrix) <- pseudobulk_ids

locked_genes <- unique(na.omit(locked$human_gene))
locked_genes <- intersect(locked_genes, rownames(count_matrix))
disease_metadata <- metadata[condition %in% c("BPD", "BPD+PH")]

expression_results <- list()
normalization_rows <- list()
result_index <- 0L
normalization_index <- 0L

for (subtype_name in unique(disease_metadata$subtype)) {
  subtype_metadata <- disease_metadata[
    subtype == subtype_name
  ][match(
    c("BPD 1", "BPD 2", "BPD+PH 1", "BPD+PH 2"),
    sample
  )]
  if (anyNA(subtype_metadata$sample)) {
    stop("Incomplete disease donor set for subtype: ", subtype_name)
  }
  subtype_counts <- count_matrix[
    ,
    subtype_metadata$pseudobulk_id,
    drop = FALSE
  ]
  dge <- DGEList(
    counts = subtype_counts,
    group = factor(
      subtype_metadata$condition,
      levels = c("BPD", "BPD+PH")
    )
  )
  dge <- calcNormFactors(dge, method = "TMM")
  log_cpm_tmm <- cpm(
    dge,
    log = TRUE,
    prior.count = 0.5,
    normalized.lib.sizes = TRUE
  )
  log_cpm_library <- cpm(
    dge,
    log = TRUE,
    prior.count = 0.5,
    normalized.lib.sizes = FALSE
  )
  cpm_tmm <- cpm(
    dge,
    log = FALSE,
    normalized.lib.sizes = TRUE
  )
  cpm_library <- cpm(
    dge,
    log = FALSE,
    normalized.lib.sizes = FALSE
  )

  for (sample_index in seq_len(nrow(subtype_metadata))) {
    sample_name <- subtype_metadata$sample[sample_index]
    result_index <- result_index + 1L
    expression_results[[result_index]] <- data.table(
      sample = sample_name,
      condition = subtype_metadata$condition[sample_index],
      subtype = subtype_name,
      cells = subtype_metadata$cells[sample_index],
      library_size = subtype_metadata$library_size[sample_index],
      gene = locked_genes,
      raw_count = subtype_counts[locked_genes, sample_index],
      detected = subtype_counts[locked_genes, sample_index] > 0,
      log2CPM_TMM = log_cpm_tmm[locked_genes, sample_index],
      log2CPM_library = log_cpm_library[locked_genes, sample_index],
      CPM_TMM = cpm_tmm[locked_genes, sample_index],
      CPM_library = cpm_library[locked_genes, sample_index]
    )
    normalization_index <- normalization_index + 1L
    normalization_rows[[normalization_index]] <- data.table(
      sample = sample_name,
      condition = subtype_metadata$condition[sample_index],
      subtype = subtype_name,
      cells = subtype_metadata$cells[sample_index],
      library_size = dge$samples$lib.size[sample_index],
      TMM_norm_factor = dge$samples$norm.factors[sample_index],
      effective_library_size = (
        dge$samples$lib.size[sample_index]
        * dge$samples$norm.factors[sample_index]
      ),
      median_metadata_nCount_RNA = (
        subtype_metadata$median_metadata_nCount_RNA[sample_index]
      ),
      median_metadata_nFeature_RNA = (
        subtype_metadata$median_metadata_nFeature_RNA[sample_index]
      ),
      median_metadata_percent_mt = (
        subtype_metadata$median_metadata_percent_mt[sample_index]
      )
    )
  }
}

expression_table <- rbindlist(expression_results)
normalization_table <- rbindlist(normalization_rows)
setorder(expression_table, subtype, sample, gene)
setorder(normalization_table, subtype, sample)

fwrite(
  expression_table,
  file.path(
    result_dir,
    "GSE275938_subtype_locked_gene_expression.tsv.gz"
  ),
  sep = "\t",
  quote = FALSE
)
fwrite(
  normalization_table,
  file.path(
    result_dir,
    "GSE275938_subtype_TMM_normalization_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(
    result_dir,
    "GSE275938_subtype_TMM_sessionInfo.txt"
  )
)

cat("Locked genes:", length(locked_genes), "\n")
cat("Expression rows:", nrow(expression_table), "\n")
print(normalization_table)
