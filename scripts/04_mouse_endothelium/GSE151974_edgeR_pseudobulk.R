#!/usr/bin/env Rscript

# Fresh animal-level endothelial pseudobulk analysis of the deposited raw UMI
# matrix. No cell is treated as an independent biological replicate.

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
  library(edgeR)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: Rscript GSE151974_edgeR_pseudobulk.R ",
    "<pseudobulk_counts.tsv.gz> <pseudobulk_metadata.tsv> <output_directory>"
  )
}

counts_path <- normalizePath(args[[1]], mustWork = TRUE)
metadata_path <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw_counts <- read_tsv(counts_path, show_col_types = FALSE)
metadata <- read_tsv(metadata_path, show_col_types = FALSE)

gene <- as.character(raw_counts[[1]])
count_matrix <- as.matrix(raw_counts[, -1])
rownames(count_matrix) <- gene
storage.mode(count_matrix) <- "integer"
count_matrix <- rowsum(count_matrix, group = rownames(count_matrix), reorder = FALSE)

if (any(count_matrix < 0L)) {
  stop("Negative counts detected.")
}
if (anyDuplicated(metadata$pseudobulk_id)) {
  stop("Duplicated pseudobulk identifiers detected.")
}
if (!setequal(colnames(count_matrix), metadata$pseudobulk_id)) {
  stop("Count columns and metadata pseudobulk identifiers do not match.")
}
metadata <- metadata[match(colnames(count_matrix), metadata$pseudobulk_id), ]
metadata$Age <- factor(metadata$Age, levels = c("P3", "P7", "P14"))
metadata$Oxygen <- factor(
  metadata$Oxygen,
  levels = c("Normoxia", "Hyperoxia")
)
metadata$CellType <- factor(
  metadata$CellType,
  levels = c("Cap", "Cap-a", "Art", "Vein", "Lymph")
)

model_specs <- list(
  source_matched_oxygen_only_all_ages = list(
    subset = function(frame) rep(TRUE, nrow(frame)),
    formula = ~ Oxygen,
    contrast_type = "oxygen_coefficient"
  ),
  age_adjusted_all_ages = list(
    subset = function(frame) rep(TRUE, nrow(frame)),
    formula = ~ Age + Oxygen,
    contrast_type = "oxygen_coefficient"
  ),
  age_by_oxygen_P14_contrast = list(
    subset = function(frame) rep(TRUE, nrow(frame)),
    formula = ~ Age * Oxygen,
    contrast_type = "p14_interaction"
  ),
  P14_oxygen_only = list(
    subset = function(frame) frame$Age == "P14",
    formula = ~ Oxygen,
    contrast_type = "oxygen_coefficient"
  ),
  age_adjusted_min10_cells = list(
    subset = function(frame) frame$n_cells >= 10,
    formula = ~ Age + Oxygen,
    contrast_type = "oxygen_coefficient"
  ),
  age_by_oxygen_P14_contrast_min10_cells = list(
    subset = function(frame) frame$n_cells >= 10,
    formula = ~ Age * Oxygen,
    contrast_type = "p14_interaction"
  ),
  P14_oxygen_only_min10_cells = list(
    subset = function(frame) frame$Age == "P14" & frame$n_cells >= 10,
    formula = ~ Oxygen,
    contrast_type = "oxygen_coefficient"
  ),
  P14_oxygen_only_min20_cells = list(
    subset = function(frame) frame$Age == "P14" & frame$n_cells >= 20,
    formula = ~ Oxygen,
    contrast_type = "oxygen_coefficient"
  )
)

all_results <- list()
summary_rows <- list()
sample_qc_rows <- list()
error_rows <- list()

for (cell_type in levels(metadata$CellType)) {
  cell_metadata <- metadata[metadata$CellType == cell_type, , drop = FALSE]
  cell_counts <- count_matrix[, cell_metadata$pseudobulk_id, drop = FALSE]

  for (model_name in names(model_specs)) {
    specification <- model_specs[[model_name]]
    keep_samples <- specification$subset(cell_metadata)
    model_metadata <- droplevels(
      cell_metadata[keep_samples, , drop = FALSE]
    )
    model_counts <- cell_counts[
      ,
      model_metadata$pseudobulk_id,
      drop = FALSE
    ]

    n_normoxia <- sum(model_metadata$Oxygen == "Normoxia")
    n_hyperoxia <- sum(model_metadata$Oxygen == "Hyperoxia")
    if (n_normoxia < 3L || n_hyperoxia < 3L) {
      error_rows[[length(error_rows) + 1L]] <- data.frame(
        CellType = cell_type,
        model = model_name,
        error = "Fewer than three pseudobulks in an oxygen group"
      )
      next
    }
    if (
      specification$contrast_type == "p14_interaction" &&
        (
          sum(
            model_metadata$Age == "P14" &
              model_metadata$Oxygen == "Normoxia"
          ) < 3L ||
            sum(
              model_metadata$Age == "P14" &
                model_metadata$Oxygen == "Hyperoxia"
            ) < 3L
        )
    ) {
      error_rows[[length(error_rows) + 1L]] <- data.frame(
        CellType = cell_type,
        model = model_name,
        error = "Fewer than three P14 pseudobulks in an oxygen group"
      )
      next
    }

    design <- model.matrix(specification$formula, data = model_metadata)
    if (qr(design)$rank < ncol(design)) {
      error_rows[[length(error_rows) + 1L]] <- data.frame(
        CellType = cell_type,
        model = model_name,
        error = "Design matrix is not full rank"
      )
      next
    }
    oxygen_coefficient <- grep("^OxygenHyperoxia$", colnames(design))
    if (length(oxygen_coefficient) != 1L) {
      stop("Could not identify OxygenHyperoxia coefficient.")
    }
    contrast <- rep(0, ncol(design))
    names(contrast) <- colnames(design)
    contrast[oxygen_coefficient] <- 1
    if (specification$contrast_type == "p14_interaction") {
      interaction_coefficient <- grep(
        "^(AgeP14:OxygenHyperoxia|OxygenHyperoxia:AgeP14)$",
        colnames(design)
      )
      if (length(interaction_coefficient) != 1L) {
        stop("Could not identify the P14-by-hyperoxia interaction coefficient.")
      }
      contrast[interaction_coefficient] <- 1
    }

    fit_result <- tryCatch(
      {
        y <- DGEList(counts = model_counts)
        keep_genes <- filterByExpr(y, design = design)
        y <- y[keep_genes, , keep.lib.sizes = FALSE]
        y <- calcNormFactors(y, method = "TMM")
        y <- estimateDisp(y, design, robust = TRUE)
        fit <- glmQLFit(y, design, robust = TRUE)
        test <- glmQLFTest(fit, contrast = contrast)
        ql_table <- topTags(test, n = Inf, sort.by = "PValue")$table
        ql_table$gene <- rownames(ql_table)
        treat <- glmTreat(fit, contrast = contrast, lfc = 1)
        treat_table <- topTags(treat, n = Inf, sort.by = "none")$table
        treat_table$gene <- rownames(treat_table)
        ql_table <- ql_table %>%
          left_join(
            treat_table %>%
              transmute(
                gene,
                PValue_treat_lfc1 = PValue,
                FDR_treat_lfc1 = FDR
              ),
            by = "gene"
          )
        ql_table
      },
      error = function(error) error
    )

    if (inherits(fit_result, "error")) {
      error_rows[[length(error_rows) + 1L]] <- data.frame(
        CellType = cell_type,
        model = model_name,
        error = conditionMessage(fit_result)
      )
      next
    }

    fit_result$CellType <- cell_type
    fit_result$model <- model_name
    fit_result$n_pseudobulks <- nrow(model_metadata)
    fit_result$n_normoxia <- n_normoxia
    fit_result$n_hyperoxia <- n_hyperoxia
    fit_result$minimum_cells_per_pseudobulk <- min(model_metadata$n_cells)
    fit_result$median_cells_per_pseudobulk <- median(model_metadata$n_cells)
    all_results[[length(all_results) + 1L]] <- fit_result[
      ,
      c(
        "gene",
        "CellType",
        "model",
        "logFC",
        "logCPM",
        "F",
        "PValue",
        "FDR",
        "PValue_treat_lfc1",
        "FDR_treat_lfc1",
        "n_pseudobulks",
        "n_normoxia",
        "n_hyperoxia",
        "minimum_cells_per_pseudobulk",
        "median_cells_per_pseudobulk"
      )
    ]

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      CellType = cell_type,
      model = model_name,
      genes_tested = nrow(fit_result),
      genes_fdr05 = sum(fit_result$FDR < 0.05),
      genes_fdr05_abs_log2fc1 = sum(
        fit_result$FDR < 0.05 & abs(fit_result$logFC) >= 1
      ),
      genes_treat_lfc1_fdr05 = sum(
        fit_result$FDR_treat_lfc1 < 0.05,
        na.rm = TRUE
      ),
      n_pseudobulks = nrow(model_metadata),
      n_normoxia = n_normoxia,
      n_hyperoxia = n_hyperoxia,
      minimum_cells_per_pseudobulk = min(model_metadata$n_cells),
      median_cells_per_pseudobulk = median(model_metadata$n_cells)
    )

    sample_qc_rows[[length(sample_qc_rows) + 1L]] <- data.frame(
      CellType = cell_type,
      model = model_name,
      pseudobulk_id = model_metadata$pseudobulk_id,
      animal_id = model_metadata$animal_id,
      Age = model_metadata$Age,
      Oxygen = model_metadata$Oxygen,
      n_cells = model_metadata$n_cells,
      library_size = colSums(model_counts),
      stringsAsFactors = FALSE
    )
  }
}

result_table <- bind_rows(all_results)
summary_table <- bind_rows(summary_rows)
sample_qc <- bind_rows(sample_qc_rows)
error_table <- bind_rows(error_rows)

write_tsv(
  result_table,
  file.path(output_dir, "GSE151974_raw_pseudobulk_edgeR_all_results.tsv.gz"),
  na = "NA"
)
write_tsv(
  summary_table,
  file.path(output_dir, "GSE151974_raw_pseudobulk_edgeR_summary.tsv"),
  na = "NA"
)
write_tsv(
  sample_qc,
  file.path(output_dir, "GSE151974_raw_pseudobulk_edgeR_sample_QC.tsv"),
  na = "NA"
)
write_tsv(
  error_table,
  file.path(output_dir, "GSE151974_raw_pseudobulk_edgeR_errors.tsv"),
  na = "NA"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "GSE151974_raw_pseudobulk_edgeR_sessionInfo.txt")
)

print(summary_table)
if (nrow(error_table)) {
  print(error_table)
}
