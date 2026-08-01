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
  library(fgsea)
})

set.seed(20260725)

exact_dir <- file.path(project_root, "07_results", "mouse_endothelium_exact")
output_dir <- file.path(project_root, "07_results", "pathway_exact")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

msigdb_dir <- file.path(
  project_root, "external_data", "msigdb", "2026.1.Mm"
)
collections <- c(
  Hallmark = "mh.all.v2026.1.Mm.symbols.gmt",
  Reactome = "m2.cp.reactome.v2026.1.Mm.symbols.gmt",
  GO_BP = "m5.go.bp.v2026.1.Mm.symbols.gmt"
)

message("Loading fixed Mouse MSigDB 2026.1 collections")
pathway_collections <- lapply(
  file.path(msigdb_dir, collections),
  fgsea::gmtPathways
)
names(pathway_collections) <- names(collections)

bulk <- fread(
  file.path(exact_dir, "GSE216046_DESeq2_all_genes.tsv"),
  na.strings = c("", "NA")
)
bulk <- bulk[
  is.finite(stat) & !is.na(gene) & gene != "",
  .(gene, rank_statistic = stat)
]
bulk <- bulk[!duplicated(gene)]

edge <- fread(
  file.path(
    exact_dir,
    "GSE151974_raw_pseudobulk_edgeR_all_results.tsv.gz"
  ),
  na.strings = c("", "NA")
)
edge <- edge[
  model %in% c(
    "age_by_oxygen_P14_contrast",
    "P14_oxygen_only",
    "age_adjusted_all_ages"
  ) &
    is.finite(F) & is.finite(logFC) & !is.na(gene) & gene != ""
]
edge[, rank_statistic := sign(logFC) * sqrt(F)]

rank_sets <- list(
  list(
    dataset = "GSE216046",
    cell_type = "bulk_lung",
    model = "DESeq2_hyperoxia_vs_air",
    ranks = bulk[, setNames(rank_statistic, gene)]
  )
)

for (model_name in c(
  "age_by_oxygen_P14_contrast",
  "P14_oxygen_only",
  "age_adjusted_all_ages"
)) {
  for (cell_type_name in sort(unique(edge[model == model_name, CellType]))) {
    subset_table <- edge[
      model == model_name & CellType == cell_type_name
    ]
    subset_table <- subset_table[!duplicated(gene)]
    rank_sets[[length(rank_sets) + 1L]] <- list(
      dataset = "GSE151974",
      cell_type = cell_type_name,
      model = model_name,
      ranks = subset_table[, setNames(rank_statistic, gene)]
    )
  }
}

run_one <- function(rank_record, collection_name, pathway_list) {
  ranks <- sort(rank_record$ranks, decreasing = TRUE)
  ranks <- ranks[is.finite(ranks) & !is.na(names(ranks)) & names(ranks) != ""]
  message(
    sprintf(
      "fgsea: %s | %s | %s | %s | %d genes",
      rank_record$dataset,
      rank_record$cell_type,
      rank_record$model,
      collection_name,
      length(ranks)
    )
  )
  result <- fgseaMultilevel(
    pathways = pathway_list,
    stats = ranks,
    minSize = 15,
    maxSize = 500,
    eps = 0,
    scoreType = "std",
    nPermSimple = 20000,
    nproc = 1
  )
  result[, `:=`(
    dataset = rank_record$dataset,
    cell_type = rank_record$cell_type,
    model = rank_record$model,
    collection = collection_name,
    n_ranked_genes = length(ranks)
  )]
  setcolorder(
    result,
    c(
      "dataset", "cell_type", "model", "collection",
      "pathway", "pval", "padj", "log2err", "ES", "NES", "size",
      "n_ranked_genes", "leadingEdge"
    )
  )
  result
}

results <- list()
result_index <- 0L
for (rank_record in rank_sets) {
  for (collection_name in names(pathway_collections)) {
    result_index <- result_index + 1L
    results[[result_index]] <- run_one(
      rank_record,
      collection_name,
      pathway_collections[[collection_name]]
    )
  }
}
all_results <- rbindlist(results, use.names = TRUE, fill = TRUE)

all_results_to_write <- copy(all_results)
all_results_to_write[
  ,
  leadingEdge := vapply(leadingEdge, paste, collapse = ";", character(1))
]
fwrite(
  all_results_to_write,
  file.path(output_dir, "mouse_fgsea_all_results.tsv.gz"),
  sep = "\t",
  quote = FALSE
)

summary_table <- all_results[
  ,
  .(
    pathways_tested = .N,
    pathways_fdr05 = sum(padj < 0.05, na.rm = TRUE),
    pathways_up_fdr05 = sum(padj < 0.05 & NES > 0, na.rm = TRUE),
    pathways_down_fdr05 = sum(padj < 0.05 & NES < 0, na.rm = TRUE)
  ),
  by = .(dataset, cell_type, model, collection)
]
setorder(summary_table, dataset, model, cell_type, collection)
fwrite(
  summary_table,
  file.path(output_dir, "mouse_fgsea_summary.tsv"),
  sep = "\t",
  quote = FALSE
)

bulk_results <- all_results[
  dataset == "GSE216046",
  .(
    collection,
    pathway,
    bulk_padj = padj,
    bulk_NES = NES,
    bulk_size = size,
    bulk_leading_edge = leadingEdge
  )
]

edge_results <- all_results[
  dataset == "GSE151974",
  .(
    collection,
    pathway,
    cell_type,
    model,
    endothelial_padj = padj,
    endothelial_NES = NES,
    endothelial_size = size,
    endothelial_leading_edge = leadingEdge
  )
]

replication <- merge(
  edge_results,
  bulk_results,
  by = c("collection", "pathway"),
  all.x = TRUE
)
replication[
  ,
  `:=`(
    both_fdr05 = !is.na(bulk_padj) & bulk_padj < 0.05 &
      endothelial_padj < 0.05,
    same_direction = !is.na(bulk_NES) &
      sign(bulk_NES) == sign(endothelial_NES)
  )
]
replication[, replicated := both_fdr05 & same_direction]
replication[, direction := fifelse(bulk_NES > 0, "up", "down")]
replication[
  ,
  leading_edge_overlap := Map(
    function(x, y) intersect(x, y),
    bulk_leading_edge,
    endothelial_leading_edge
  )
]
replication[, leading_edge_overlap_n := lengths(leading_edge_overlap)]

replication_to_write <- copy(replication)
for (column_name in c(
  "bulk_leading_edge",
  "endothelial_leading_edge",
  "leading_edge_overlap"
)) {
  replication_to_write[
    ,
    (column_name) := vapply(
      get(column_name), paste, collapse = ";", character(1)
    )
  ]
}
fwrite(
  replication_to_write,
  file.path(output_dir, "mouse_pathway_replication_all.tsv.gz"),
  sep = "\t",
  quote = FALSE
)

primary_replication <- replication[
  model == "age_by_oxygen_P14_contrast" & replicated
]
setorder(
  primary_replication,
  collection,
  pathway,
  -leading_edge_overlap_n,
  cell_type
)
primary_to_write <- copy(primary_replication)
for (column_name in c(
  "bulk_leading_edge",
  "endothelial_leading_edge",
  "leading_edge_overlap"
)) {
  primary_to_write[
    ,
    (column_name) := vapply(
      get(column_name), paste, collapse = ";", character(1)
    )
  ]
}
fwrite(
  primary_to_write,
  file.path(output_dir, "mouse_pathway_replication_primary.tsv"),
  sep = "\t",
  quote = FALSE
)

p14_support <- replication[
  model == "P14_oxygen_only",
  .(
    collection,
    pathway,
    cell_type,
    p14_replicated = replicated,
    p14_padj = endothelial_padj,
    p14_NES = endothelial_NES
  )
]

age_adjusted_support <- replication[
  model == "age_adjusted_all_ages",
  .(
    collection,
    pathway,
    cell_type,
    age_adjusted_replicated = replicated,
    age_adjusted_padj = endothelial_padj,
    age_adjusted_NES = endothelial_NES
  )
]

pathway_breadth <- primary_replication[
  ,
  .(
    direction = first(direction),
    primary_subtypes_replicated = uniqueN(cell_type),
    primary_subtypes = paste(sort(unique(cell_type)), collapse = ";"),
    bulk_NES = first(bulk_NES),
    bulk_padj = first(bulk_padj),
    median_endothelial_NES = median(endothelial_NES, na.rm = TRUE),
    minimum_leading_edge_overlap = as.numeric(
      min(leading_edge_overlap_n)
    ),
    median_leading_edge_overlap = as.numeric(
      median(leading_edge_overlap_n)
    )
  ),
  by = .(collection, pathway)
]

p14_breadth <- p14_support[
  p14_replicated == TRUE,
  .(
    p14_subtypes_replicated = uniqueN(cell_type),
    p14_subtypes = paste(sort(unique(cell_type)), collapse = ";")
  ),
  by = .(collection, pathway)
]
pathway_breadth <- merge(
  pathway_breadth,
  p14_breadth,
  by = c("collection", "pathway"),
  all.x = TRUE
)
age_adjusted_breadth <- age_adjusted_support[
  age_adjusted_replicated == TRUE,
  .(
    age_adjusted_subtypes_replicated = uniqueN(cell_type),
    age_adjusted_subtypes = paste(sort(unique(cell_type)), collapse = ";")
  ),
  by = .(collection, pathway)
]
pathway_breadth <- merge(
  pathway_breadth,
  age_adjusted_breadth,
  by = c("collection", "pathway"),
  all.x = TRUE
)
pathway_breadth[
  is.na(p14_subtypes_replicated),
  `:=`(p14_subtypes_replicated = 0L, p14_subtypes = "")
]
pathway_breadth[
  is.na(age_adjusted_subtypes_replicated),
  `:=`(
    age_adjusted_subtypes_replicated = 0L,
    age_adjusted_subtypes = ""
  )
]
setorder(
  pathway_breadth,
  -primary_subtypes_replicated,
  -p14_subtypes_replicated,
  -age_adjusted_subtypes_replicated,
  bulk_padj
)
fwrite(
  pathway_breadth,
  file.path(output_dir, "mouse_pathway_replication_breadth.tsv"),
  sep = "\t",
  quote = FALSE
)

replication_summary <- replication[
  ,
  .(
    pathways_compared = .N,
    pathways_both_fdr05 = sum(both_fdr05, na.rm = TRUE),
    pathways_replicated_same_direction = sum(replicated, na.rm = TRUE),
    pathways_replicated_up = sum(replicated & direction == "up", na.rm = TRUE),
    pathways_replicated_down = sum(replicated & direction == "down", na.rm = TRUE)
  ),
  by = .(model, cell_type, collection)
]
setorder(replication_summary, model, cell_type, collection)
fwrite(
  replication_summary,
  file.path(output_dir, "mouse_pathway_replication_summary.tsv"),
  sep = "\t",
  quote = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "mouse_pathway_fgsea_sessionInfo.txt")
)

message("Pathway analysis complete")
print(summary_table)
print(head(pathway_breadth, 25))
