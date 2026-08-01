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
.libPaths(c(file.path(project_root, ".R_library"), .libPaths()))

suppressPackageStartupMessages(library(data.table))

result_root <- file.path(project_root, "07_results")
figure_dir <- file.path(project_root, "08_figures", "manuscript_main")
plot_data_dir <- file.path(project_root, "09_plotting_data", "manuscript_main")

checks <- list()
add_check <- function(check, expected, observed, pass) {
  checks[[length(checks) + 1L]] <<- data.table(
    check = check,
    expected = as.character(expected),
    observed = as.character(observed),
    pass = isTRUE(pass)
  )
}

bulk <- fread(file.path(
  result_root,
  "mouse_endothelium_exact",
  "GSE216046_DESeq2_all_genes.tsv"
))
bulk_threshold <- bulk[
  !is.na(padj)
    & padj < 0.05
    & abs(log2FoldChange_mle) >= 1
]
add_check("Figure 2 genes after filter", 14561, nrow(bulk), nrow(bulk) == 14561)
add_check(
  "Figure 2 threshold-positive genes",
  1029,
  nrow(bulk_threshold),
  nrow(bulk_threshold) == 1029
)
bulk_lfc_threshold <- bulk[
  !is.na(padj_lfcThreshold1) & padj_lfcThreshold1 < 0.05
]
add_check(
  "Figure 2 DESeq2 lfcThreshold=1 genes",
  341,
  nrow(bulk_lfc_threshold),
  nrow(bulk_lfc_threshold) == 341
)

pca <- fread(file.path(
  result_root,
  "mouse_endothelium_exact",
  "GSE216046_DESeq2_vst_PCA.tsv"
))
pca_counts <- pca[, .N, by = oxygen]
add_check(
  "Figure 2 PCA samples per group",
  "Air=4;Hyperoxia=4",
  paste0(pca_counts$oxygen, "=", pca_counts$N, collapse = ";"),
  identical(pca_counts$N, c(4L, 4L))
)

convergent <- fread(file.path(
  result_root,
  "mouse_endothelium_exact",
  "mouse_endothelial_convergent_genes.tsv"
))[model == "age_by_oxygen_P14_contrast"]
locked_gene_count <- uniqueN(convergent[direction_concordant == TRUE, toupper(gene)])
add_check(
  "Figure 1 sequentially defined direction-concordant genes",
  198,
  locked_gene_count,
  locked_gene_count == 198
)

pathway_primary <- fread(file.path(
  result_root,
  "pathway_exact",
  "mouse_pathway_replication_primary.tsv"
))
replicated_pathway_count <- uniqueN(pathway_primary$pathway)
add_check(
  "Figure 1 unique primary replicated pathways",
  162,
  replicated_pathway_count,
  replicated_pathway_count == 162
)

composition_tests <- fread(file.path(
  result_root,
  "mouse_endothelium_exact",
  "GSE151974_composition_logit_exact_permutation.tsv"
))
p14_gcap <- composition_tests[
  Age == "P14"
    & metric == "Cap_fraction_of_endothelium",
  100 * raw_fraction_difference_hyperoxia_minus_normoxia
]
p14_acap <- composition_tests[
  Age == "P14"
    & metric == "Cap-a_fraction_of_endothelium",
  100 * raw_fraction_difference_hyperoxia_minus_normoxia
]
add_check(
  "Figure 3 P14 gCap-like difference (percentage points)",
  "-30.705 (+/-0.001)",
  sprintf("%.3f", p14_gcap),
  abs(p14_gcap + 30.7048163) < 0.001
)
add_check(
  "Figure 3 P14 aCap-like difference (percentage points)",
  "29.358 (+/-0.001)",
  sprintf("%.3f", p14_acap),
  abs(p14_acap - 29.3579428) < 0.001
)

replication <- fread(file.path(
  result_root,
  "mouse_endothelium_exact",
  "mouse_endothelial_cross_dataset_replication_summary.tsv"
))[model == "age_by_oxygen_P14_contrast"]
cap <- replication[CellType == "Cap"]
cap_a <- replication[CellType == "Cap-a"]
add_check(
  "Figure 4 Cap replication",
  "174/179 concordant; rho=0.653",
  sprintf(
    "%d/%d concordant; rho=%.3f",
    cap$direction_concordant_overlap_genes,
    cap$significant_overlap_genes,
    cap$spearman_all_shared_gene_effects
  ),
  cap$direction_concordant_overlap_genes == 174
    && cap$significant_overlap_genes == 179
    && abs(cap$spearman_all_shared_gene_effects - 0.652870) < 0.0001
)
add_check(
  "Figure 4 Cap-a replication",
  "127/128 concordant; rho=0.598",
  sprintf(
    "%d/%d concordant; rho=%.3f",
    cap_a$direction_concordant_overlap_genes,
    cap_a$significant_overlap_genes,
    cap_a$spearman_all_shared_gene_effects
  ),
  cap_a$direction_concordant_overlap_genes == 127
    && cap_a$significant_overlap_genes == 128
    && abs(cap_a$spearman_all_shared_gene_effects - 0.597795) < 0.0001
)

human_core <- fread(file.path(
  result_root,
  "human_endothelial_subtypes",
  "GSE275938_capillary_robust_core_genes.tsv"
))
add_check(
  "Figure 6 exploratory post-projection genes",
  17,
  nrow(human_core),
  nrow(human_core) == 17
)

orthologs <- fread(file.path(
  result_root,
  "cross_species",
  "formal_one_to_one_mouse_signature.tsv"
))
ccn2_mapping <- orthologs[
  gene == "Ctgf" & human_gene == "CCN2"
]
add_check(
  "CCN2/Ctgf formal NCBI synonym recovery",
  "one-to-one; unique NCBI synonym",
  paste0(
    ccn2_mapping$orthology_cardinality,
    "; ",
    ccn2_mapping$mouse_symbol_resolution
  ),
  nrow(ccn2_mapping) == 1
    && ccn2_mapping$orthology_cardinality == "one-to-one"
    && ccn2_mapping$mouse_symbol_resolution
      == "unique match in NCBI Synonyms"
)

human_compact <- fread(file.path(
  result_root,
  "human_endothelial_subtypes",
  "GSE275938_human_subtype_compact_evidence.tsv"
))
human_capillary <- human_compact[subtype %in% c("gCap", "aCap")]
add_check(
  "Figure 5 human capillary empirical BH q-values",
  "gCap<0.01;aCap<0.01",
  paste0(
    human_capillary$subtype,
    "=",
    sprintf(
      "%.4f",
      human_capillary$empirical_q_greater_BH_within_signature
    ),
    collapse = ";"
  ),
  all(
    human_capillary$empirical_q_greater_BH_within_signature < 0.01
  )
)

clr <- fread(file.path(
  result_root,
  "mouse_endothelium_exact",
  "GSE151974_composition_CLR_exact_permutation.tsv"
))
p14_clr_p <- unique(clr[Age == "P14", global_exact_permutation_p])
add_check(
  "Figure 6 P14 CLR exact permutation",
  "0.002165",
  sprintf("%.6f", p14_clr_p),
  length(p14_clr_p) == 1
    && abs(p14_clr_p - 0.002164502) < 1e-6
)

expected_stems <- c(
  "Figure1_study_design",
  "Figure2_GSE216046_discovery",
  "Figure3_GSE151974_composition",
  "Figure4_mouse_replication",
  "Figure5_human_subtype_projection",
  "Figure6_core_and_model"
)
for (stem in expected_stems) {
  for (extension in c("png", "pdf")) {
    path <- file.path(figure_dir, paste0(stem, ".", extension))
    observed_size <- if (file.exists(path)) file.info(path)$size else 0
    add_check(
      paste0("File integrity: ", stem, ".", extension),
      "exists and >5 KB",
      paste0(observed_size, " bytes"),
      file.exists(path) && observed_size > 5000
    )
  }
}

qa <- rbindlist(checks)
fwrite(
  qa,
  file.path(figure_dir, "manuscript_main_figures_QA.tsv"),
  sep = "\t"
)

if (!all(qa$pass)) {
  print(qa[pass == FALSE])
  stop("One or more manuscript-figure QA checks failed.")
}

cat(sprintf("All %d manuscript-figure QA checks passed.\n", nrow(qa)))
