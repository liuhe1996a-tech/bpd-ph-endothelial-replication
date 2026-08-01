#!/usr/bin/env Rscript

# Exact count-model analysis planned for the manuscript.
# Requires R, DESeq2, and readr. Packages are loaded from the project-local
# library so the user's system R library remains unchanged.

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
  library(DESeq2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: Rscript GSE216046_DESeq2.R ",
    "<GSE216046_gene_count.csv.gz> <output_directory>"
  )
}

input_path <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sample_names <- c(
  "Air_1", "Air_2", "Air_3", "Air_4",
  "Hyperoxia_1", "Hyperoxia_2", "Hyperoxia_3", "Hyperoxia_4"
)

raw <- read_csv(input_path, show_col_types = FALSE)
raw$gene_name[is.na(raw$gene_name) | raw$gene_name == ""] <- raw$gene_id[
  is.na(raw$gene_name) | raw$gene_name == ""
]

count_matrix <- rowsum(
  as.matrix(raw[, sample_names]),
  group = raw$gene_name,
  reorder = FALSE
)
storage.mode(count_matrix) <- "integer"

coldata <- data.frame(
  row.names = sample_names,
  oxygen = factor(
    c(rep("Air", 4), rep("Hyperoxia", 4)),
    levels = c("Air", "Hyperoxia")
  )
)

# Independent low-count filtering; this matches the pre-specified Python
# sensitivity analysis and is fixed before model fitting.
library_sizes <- colSums(count_matrix)
cpm <- sweep(count_matrix, 2L, library_sizes, "/") * 1e6
keep <- rowSums(cpm >= 1) >= 3L
count_matrix <- count_matrix[keep, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ oxygen
)
dds <- DESeq(dds, betaPrior = FALSE, quiet = TRUE)

res_unshrunken <- results(
  dds,
  contrast = c("oxygen", "Hyperoxia", "Air"),
  alpha = 0.05,
  independentFiltering = TRUE
)
res_lfc_threshold_1 <- results(
  dds,
  contrast = c("oxygen", "Hyperoxia", "Air"),
  alpha = 0.05,
  lfcThreshold = 1,
  altHypothesis = "greaterAbs",
  independentFiltering = TRUE
)
res_shrunken <- lfcShrink(
  dds,
  coef = "oxygen_Hyperoxia_vs_Air",
  type = "normal"
)

result <- data.frame(
  gene = rownames(res_unshrunken),
  baseMean = res_unshrunken$baseMean,
  log2FoldChange_mle = res_unshrunken$log2FoldChange,
  lfcSE_mle = res_unshrunken$lfcSE,
  stat = res_unshrunken$stat,
  pvalue = res_unshrunken$pvalue,
  padj = res_unshrunken$padj,
  stat_lfcThreshold1 = res_lfc_threshold_1$stat,
  pvalue_lfcThreshold1 = res_lfc_threshold_1$pvalue,
  padj_lfcThreshold1 = res_lfc_threshold_1$padj,
  log2FoldChange_normal_shrunk = res_shrunken$log2FoldChange,
  lfcSE_normal_shrunk = res_shrunken$lfcSE,
  stringsAsFactors = FALSE
)
result <- result[order(result$padj, result$pvalue, na.last = TRUE), ]

write_tsv(
  result,
  file.path(output_dir, "GSE216046_DESeq2_all_genes.tsv"),
  na = "NA"
)
write_tsv(
  subset(
    result,
    !is.na(padj) & padj < 0.05 & abs(log2FoldChange_mle) >= 1
  ),
  file.path(output_dir, "GSE216046_DESeq2_fdr05_abs_mle_log2fc1.tsv"),
  na = "NA"
)
write_tsv(
  subset(
    result,
    !is.na(padj_lfcThreshold1) & padj_lfcThreshold1 < 0.05
  ),
  file.path(output_dir, "GSE216046_DESeq2_lfcThreshold1_fdr05.tsv"),
  na = "NA"
)

sample_qc <- data.frame(
  sample = sample_names,
  oxygen = coldata$oxygen,
  library_size = library_sizes[sample_names],
  size_factor = sizeFactors(dds)[sample_names],
  detected_genes_raw = colSums(count_matrix > 0),
  stringsAsFactors = FALSE
)
write_tsv(
  sample_qc,
  file.path(output_dir, "GSE216046_DESeq2_sample_QC.tsv"),
  na = "NA"
)

vsd <- vst(dds, blind = FALSE)
pca <- plotPCA(vsd, intgroup = "oxygen", returnData = TRUE)
pca$sample <- rownames(pca)
write_tsv(
  pca,
  file.path(output_dir, "GSE216046_DESeq2_vst_PCA.tsv"),
  na = "NA"
)

session <- capture.output(sessionInfo())
writeLines(session, file.path(output_dir, "GSE216046_DESeq2_sessionInfo.txt"))

summary_lines <- c(
  paste0("genes_after_cpm_filter\t", nrow(count_matrix)),
  paste0(
    "genes_fdr05\t",
    sum(!is.na(result$padj) & result$padj < 0.05)
  ),
  paste0(
    "genes_fdr05_abs_mle_log2fc1\t",
    sum(
      !is.na(result$padj) &
        result$padj < 0.05 &
        abs(result$log2FoldChange_mle) >= 1
    )
  ),
  paste0(
    "genes_lfcThreshold1_fdr05\t",
    sum(
      !is.na(result$padj_lfcThreshold1) &
        result$padj_lfcThreshold1 < 0.05
    )
  )
)
writeLines(
  c("metric\tvalue", summary_lines),
  file.path(output_dir, "GSE216046_DESeq2_summary.tsv")
)
