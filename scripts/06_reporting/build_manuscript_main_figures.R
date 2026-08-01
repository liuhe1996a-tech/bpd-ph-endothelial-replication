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
  library(ggplot2)
  library(cowplot)
  library(grid)
})

result_root <- file.path(project_root, "07_results")
figure_dir <- file.path(project_root, "08_figures", "manuscript_main")
plot_data_dir <- file.path(project_root, "09_plotting_data", "manuscript_main")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_data_dir, recursive = TRUE, showWarnings = FALSE)

ink <- "#27313A"
muted <- "#64717C"
grid_colour <- "#D9DEE3"
air_colour <- "#2C5F8A"
hyperoxia_colour <- "#D97932"
up_colour <- "#C95046"
down_colour <- "#3F79A8"
neutral_colour <- "#C7CDD2"
gold_colour <- "#C69C3B"
green_colour <- "#4F8C77"
purple_colour <- "#755B9C"

theme_research <- theme_minimal(base_size = 10) +
  theme(
    text = element_text(colour = ink, family = "sans"),
    axis.text = element_text(colour = ink),
    axis.title = element_text(colour = ink),
    plot.title = element_text(face = "bold", colour = ink, size = 11),
    plot.subtitle = element_text(colour = muted, size = 8.8),
    panel.grid.major = element_line(colour = grid_colour, linewidth = 0.35),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 8.5),
    legend.text = element_text(size = 8.5),
    strip.text = element_text(face = "bold", colour = ink),
    plot.margin = margin(7, 10, 7, 7)
  )

save_main_figure <- function(plot, stem, width, height, dpi = 320) {
  ggsave(
    file.path(figure_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
  ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    bg = "white",
    limitsize = FALSE
  )
}

add_figure_header <- function(plot, title, subtitle, title_size = 15) {
  header <- ggdraw() +
    draw_label(
      title,
      x = 0.01,
      y = 0.75,
      hjust = 0,
      fontface = "bold",
      size = title_size,
      colour = ink
    ) +
    draw_label(
      subtitle,
      x = 0.01,
      y = 0.16,
      hjust = 0,
      size = 9.2,
      colour = muted
    )
  plot_grid(header, plot, ncol = 1, rel_heights = c(0.075, 1))
}

clean_pathway_label <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^REACTOME_", "", x)
  x <- sub("^GOBP_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  replacements <- c(
    "Oxidative Phosphorylation" = "Oxidative phosphorylation",
    "Respiratory Electron Transport" = "Respiratory electron transport",
    "Complex I Biogenesis" = "Complex I biogenesis",
    "Nuclear Events Mediated by Nfe2l2" = "NFE2L2-mediated nuclear events",
    "Mtorc1 Signaling" = "mTORC1 signaling",
    "P53 Pathway" = "p53 pathway",
    "Reactive Oxygen Species Pathway" = "Reactive oxygen species pathway",
    "Tnfa Signaling via Nfkb" = "TNF-alpha signaling via NF-kappaB",
    "Epithelial Mesenchymal Transition" = "EMT-like program",
    "Extracellular Matrix Organization" = "Extracellular matrix organization",
    "Degradation of the Extracellular Matrix" = "Extracellular matrix degradation"
  )
  matched <- match(x, names(replacements))
  x[!is.na(matched)] <- unname(replacements[matched[!is.na(matched)]])
  x
}

canonical_gene_symbol <- function(x) {
  x <- toupper(x)
  x[x == "CTGF"] <- "CCN2"
  x
}

# -------------------------------------------------------------------------
# Shared source data
# -------------------------------------------------------------------------

bulk_de <- fread(file.path(
  result_root,
  "mouse_endothelium_exact",
  "GSE216046_DESeq2_all_genes.tsv"
))
bulk_de[, gene_upper := toupper(gene)]
bulk_de[, gene_canonical := canonical_gene_symbol(gene)]

bulk_pca <- fread(file.path(
  result_root,
  "mouse_endothelium_exact",
  "GSE216046_DESeq2_vst_PCA.tsv"
))

pseudobulk <- fread(file.path(
  result_root,
  "mouse_endothelium_exact",
  "GSE151974_raw_pseudobulk_edgeR_all_results.tsv.gz"
))
pseudobulk[, gene_upper := toupper(gene)]
pseudobulk[, gene_canonical := canonical_gene_symbol(gene)]

replication_summary <- fread(file.path(
  result_root,
  "mouse_endothelium_exact",
  "mouse_endothelial_cross_dataset_replication_summary.tsv"
))[model == "age_by_oxygen_P14_contrast"]

convergent_genes <- fread(file.path(
  result_root,
  "mouse_endothelium_exact",
  "mouse_endothelial_convergent_genes.tsv"
))[model == "age_by_oxygen_P14_contrast"]
convergent_genes[, gene_upper := toupper(gene)]
convergent_genes[, gene_canonical := canonical_gene_symbol(gene)]

composition <- fread(file.path(
  project_root,
  "03_metadata",
  "GSE151974_animal_level_cell_composition.tsv"
))

composition_tests <- fread(file.path(
  result_root,
  "mouse_endothelium_exact",
  "GSE151974_composition_logit_exact_permutation.tsv"
))

fgsea_all <- fread(file.path(
  result_root,
  "pathway_exact",
  "mouse_fgsea_all_results.tsv.gz"
))

pathway_priority <- fread(file.path(
  result_root,
  "pathway_exact",
  "mouse_pathway_priority.tsv"
))

human_core <- fread(file.path(
  result_root,
  "human_endothelial_subtypes",
  "GSE275938_capillary_robust_core_genes.tsv"
))
human_core[, gene_upper := toupper(gene)]

human_compact <- fread(file.path(
  result_root,
  "human_endothelial_subtypes",
  "GSE275938_human_subtype_compact_evidence.tsv"
))

core_genes <- human_core$gene_upper
core_labels <- human_core$gene_upper

# -------------------------------------------------------------------------
# Figure 1. Study design and evidence locking
# -------------------------------------------------------------------------

workflow <- data.table(
  stage = 1:4,
  x = 1:4,
  heading = c(
    "Discovery",
    "Independent mouse replication",
    "Mouse signature definition",
    "Human projection"
  ),
  body = c(
    "GSE216046\nPurified lung endothelial bulk RNA-seq\nAir n=4 | Hyperoxia n=4\nDESeq2 + ranked GSEA",
    "GSE151974 lung scRNA-seq\n36 animals | 5 endothelial subtypes\nAnimal x subtype pseudobulk\nedgeR: categorical Age x Oxygen\nP14 hyperoxia contrast",
    "198 direction-concordant genes\n162 replicated pathway annotations\nSequentially frozen for current projection",
    "GSE275938 lung scRNA-seq\nBPD n=2 | BPD+PH n=2\n187 NCBI one-to-one orthologs\nExploratory donor x subtype projection"
  ),
  box_fill = c("#E4EEF5", "#E4EEF5", "#F5EBD4", "#E8F0EA"),
  border = c(air_colour, air_colour, gold_colour, green_colour)
)

evidence_rules <- data.table(
  x = c(1.3, 2.5, 3.7),
  label = c(
    "Biological unit:\nsample / animal / donor",
    "Discovery: FDR<0.05 and observed |MLE log2FC|>=1\nReplication: FDR<0.05 and concordant direction",
    "Human result:\nexploratory directional support"
  ),
  fill = c("#F2F5F7", "#F2F5F7", "#FFF2EB")
)

f1 <- ggplot() +
  geom_rect(
    data = workflow,
    aes(
      xmin = x - 0.42, xmax = x + 0.42,
      ymin = 0.38, ymax = 0.90,
      fill = box_fill, colour = border
    ),
    linewidth = 1.1,
    show.legend = FALSE
  ) +
  scale_fill_identity() +
  scale_colour_identity() +
  geom_text(
    data = workflow,
    aes(x = x, y = 0.83, label = heading),
    fontface = "bold",
    colour = ink,
    size = 4.0,
    lineheight = 0.95
  ) +
  geom_text(
    data = workflow,
    aes(x = x, y = 0.62, label = body),
    colour = ink,
    size = 3.25,
    lineheight = 1.12
  ) +
  geom_segment(
    data = data.table(x = c(1.43, 2.43, 3.43), xend = c(1.57, 2.57, 3.57)),
    aes(x = x, xend = xend, y = 0.64, yend = 0.64),
    linewidth = 1.0,
    colour = muted,
    arrow = arrow(length = unit(0.15, "inches"), type = "closed")
  ) +
  geom_rect(
    data = evidence_rules,
    aes(
      xmin = x - 0.48, xmax = x + 0.48,
      ymin = 0.08, ymax = 0.27,
      fill = fill
    ),
    colour = "#AEB7BE",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  geom_text(
    data = evidence_rules,
    aes(x = x, y = 0.175, label = label),
    colour = ink,
    size = 3.25,
    lineheight = 1.02
  ) +
  annotate(
    "text",
    x = 2.5,
    y = 0.985,
    label = "Sequential design reduces cell-level pseudoreplication and post-hoc analytical flexibility",
    colour = muted,
    size = 3.6
  ) +
  coord_cartesian(xlim = c(0.45, 4.55), ylim = c(0, 1.05), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(12, 20, 8, 20))

f1_final <- add_figure_header(
  f1,
  "Figure 1. Sequential cross-dataset mouse-to-human analysis design",
  "Two independent mouse datasets define and freeze a P14-focused signature before an exploratory human BPD-PH projection."
)
save_main_figure(f1_final, "Figure1_study_design", 15, 7.8)
fwrite(workflow, file.path(plot_data_dir, "Figure1_workflow.tsv"), sep = "\t")
fwrite(
  evidence_rules,
  file.path(plot_data_dir, "Figure1_evidence_rules.tsv"),
  sep = "\t"
)

# -------------------------------------------------------------------------
# Figure 2. Purified endothelial discovery dataset
# -------------------------------------------------------------------------

bulk_pca[, oxygen := factor(oxygen, levels = c("Air", "Hyperoxia"))]
p2a <- ggplot(
  bulk_pca,
  aes(x = PC1, y = PC2, colour = oxygen, shape = oxygen)
) +
  geom_hline(yintercept = 0, colour = grid_colour, linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = grid_colour, linewidth = 0.35) +
  geom_point(size = 3.8, stroke = 0.9) +
  geom_text(
    aes(label = sample),
    nudge_y = 0.85,
    size = 2.7,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c("Air" = air_colour, "Hyperoxia" = hyperoxia_colour)) +
  scale_shape_manual(values = c("Air" = 16, "Hyperoxia" = 17)) +
  labs(
    title = "A. Variance-stabilized PCA",
    subtitle = "PC1=67.67%; PC2=8.51%",
    x = "PC1",
    y = "PC2",
    colour = NULL,
    shape = NULL
  ) +
  theme_research

bulk_de[, de_class := fifelse(
  !is.na(padj) & padj < 0.05 & log2FoldChange_mle >= 1,
  "Up",
  fifelse(
    !is.na(padj) & padj < 0.05 & log2FoldChange_mle <= -1,
    "Down",
    "Not threshold-positive"
  )
)]
bulk_de[, neglog10_fdr := -log10(pmax(padj, 1e-60))]
bulk_de[, is_core := FALSE]
bulk_de[, label_gene := NA_character_]

p2b <- ggplot(
  bulk_de[!is.na(padj)],
  aes(x = log2FoldChange_mle, y = neglog10_fdr)
) +
  geom_point(
    aes(colour = de_class),
    size = 0.75,
    alpha = 0.48
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = muted, linewidth = 0.45) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = muted, linewidth = 0.45) +
  scale_colour_manual(values = c(
    "Up" = up_colour,
    "Down" = down_colour,
    "Not threshold-positive" = neutral_colour
  )) +
  labs(
    title = "B. DESeq2 differential expression",
    subtitle = "1,029/14,561 genes: FDR<0.05 and |MLE log2FC|>=1; values below FDR 1e-60 are plotted at 60",
    x = "Hyperoxia vs air MLE log2 fold-change",
    y = "-log10(FDR)",
    colour = NULL
  ) +
  theme_research +
  theme(legend.position = "top")

effect_threshold_summary <- data.table(
  definition = factor(
    c(
      "FDR < 0.05",
      "FDR < 0.05 + |observed log2FC| >= 1",
      "DESeq2 lfcThreshold = 1; FDR < 0.05"
    ),
    levels = rev(c(
      "FDR < 0.05",
      "FDR < 0.05 + |observed log2FC| >= 1",
      "DESeq2 lfcThreshold = 1; FDR < 0.05"
    ))
  ),
  genes = c(5543, 1029, 341),
  test_type = c(
    "Zero-effect test",
    "Descriptive effect filter",
    "Inferential effect-threshold test"
  )
)

p2c <- ggplot(
  effect_threshold_summary,
  aes(y = definition, x = genes, fill = test_type)
) +
  geom_col(width = 0.62, alpha = 0.92) +
  geom_text(
    aes(label = format(genes, big.mark = ",")),
    hjust = -0.1,
    colour = ink,
    size = 3.0,
    fontface = "bold"
  ) +
  scale_fill_manual(values = c(
    "Zero-effect test" = "#9CB6C8",
    "Descriptive effect filter" = gold_colour,
    "Inferential effect-threshold test" = green_colour
  )) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "C. Effect-threshold sensitivity",
    subtitle = "The |observed log2FC| filter is descriptive; lfcThreshold=1 tests effects beyond one log2 unit",
    x = "Genes meeting each definition",
    y = NULL,
    fill = NULL
  ) +
  theme_research +
  theme(legend.position = "bottom")

selected_hallmark <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)
bulk_hallmark <- fgsea_all[
  dataset == "GSE216046"
    & collection == "Hallmark"
    & pathway %in% selected_hallmark,
  .(pathway, NES, padj)
]
bulk_hallmark[, label := clean_pathway_label(pathway)]
bulk_hallmark[, label := factor(
  label,
  levels = clean_pathway_label(rev(selected_hallmark))
)]
bulk_hallmark[, fdr_label := fifelse(padj < 0.001, "FDR<0.001", sprintf("FDR=%.3f", padj))]

p2d <- ggplot(bulk_hallmark, aes(x = NES, y = label)) +
  geom_vline(xintercept = 0, colour = muted, linewidth = 0.45) +
  geom_col(aes(fill = NES > 0), width = 0.68, alpha = 0.9) +
  geom_text(
    aes(label = fdr_label),
    hjust = -0.08,
    colour = ink,
    size = 2.55
  ) +
  scale_fill_manual(values = c("TRUE" = up_colour, "FALSE" = down_colour), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(
    title = "D. Ranked Hallmark enrichment",
    subtitle = "Full-gene Wald-statistic ranking; positive NES indicates hyperoxia enrichment",
    x = "Normalized enrichment score",
    y = NULL
  ) +
  theme_research +
  theme(legend.position = "none")

f2_body <- plot_grid(
  p2a,
  p2b,
  p2c,
  p2d,
  ncol = 2,
  rel_widths = c(0.94, 1.06),
  rel_heights = c(0.94, 1.06),
  align = "hv"
)
f2_final <- add_figure_header(
  f2_body,
  "Figure 2. Purified mouse lung endothelium mounts a broad hyperoxia response",
  "GSE216046 discovery analysis uses exact raw-count DESeq2 and full-ranking enrichment (Air n=4; Hyperoxia n=4)."
)
save_main_figure(f2_final, "Figure2_GSE216046_discovery", 15.5, 11.5)
bulk_pca[, `:=`(
  PC1_explained_variance_percent = 67.67,
  PC2_explained_variance_percent = 8.51
)]
fwrite(
  bulk_pca,
  file.path(plot_data_dir, "Figure2_PCA_coordinates.tsv"),
  sep = "\t"
)
fwrite(
  bulk_de[, .(
    gene,
    gene_canonical,
    baseMean,
    log2FoldChange_mle,
    lfcSE_mle,
    padj,
    de_class,
    is_core
  )],
  file.path(plot_data_dir, "Figure2_DESeq2_volcano.tsv"),
  sep = "\t"
)
fwrite(
  effect_threshold_summary,
  file.path(plot_data_dir, "Figure2_effect_threshold_sensitivity.tsv"),
  sep = "\t"
)
fwrite(bulk_hallmark, file.path(plot_data_dir, "Figure2_selected_hallmark.tsv"), sep = "\t")

# -------------------------------------------------------------------------
# Figure 3. Animal-level endothelial composition
# -------------------------------------------------------------------------

composition[, Age := factor(Age, levels = c("P3", "P7", "P14"))]
composition[, Oxygen := factor(Oxygen, levels = c("Normoxia", "Hyperoxia"))]

composition_metrics <- c(
  "capillary_general_like_fraction_of_endothelium",
  "capillary_aerocyte_like_fraction_of_endothelium",
  "arterial_fraction_of_endothelium",
  "venous_fraction_of_endothelium",
  "lymphatic_fraction_of_endothelium"
)
composition_labels <- c(
  "capillary_general_like_fraction_of_endothelium" = "gCap-like",
  "capillary_aerocyte_like_fraction_of_endothelium" = "aCap-like",
  "arterial_fraction_of_endothelium" = "Arterial",
  "venous_fraction_of_endothelium" = "Venous",
  "lymphatic_fraction_of_endothelium" = "Lymphatic",
  "endothelial_fraction_all_cells" = "Endothelial / all cells"
)
composition_palette <- c(
  "gCap-like" = "#4D79A7",
  "aCap-like" = "#E07A5F",
  "Arterial" = "#76A56A",
  "Venous" = "#9C76A5",
  "Lymphatic" = "#D2AC52"
)

composition_long <- melt(
  composition,
  id.vars = c("animal_id", "Age", "Oxygen"),
  measure.vars = composition_metrics,
  variable.name = "metric",
  value.name = "fraction"
)
composition_long[, subtype := composition_labels[as.character(metric)]]
composition_long[, subtype := factor(
  subtype,
  levels = c("gCap-like", "aCap-like", "Arterial", "Venous", "Lymphatic")
)]
composition_means <- composition_long[, .(
  mean_fraction = mean(fraction),
  sd_fraction = sd(fraction),
  n_animals = .N
), by = .(Age, Oxygen, subtype)]

p3a <- ggplot(
  composition_means,
  aes(x = Oxygen, y = mean_fraction, fill = subtype)
) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.35) +
  facet_grid(. ~ Age) +
  scale_fill_manual(values = composition_palette) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    title = "A. Mean endothelial-subtype composition",
    subtitle = "Five subtypes sum to total recovered endothelium; n=6 animals per condition and age",
    x = NULL,
    y = "Mean fraction of endothelial cells",
    fill = NULL
  ) +
  theme_research +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

capillary_long <- composition_long[
  metric %in% c(
    "capillary_general_like_fraction_of_endothelium",
    "capillary_aerocyte_like_fraction_of_endothelium"
  )
]
capillary_long[, metric_label := factor(
  composition_labels[as.character(metric)],
  levels = c("gCap-like", "aCap-like")
)]

p3b <- ggplot(
  capillary_long,
  aes(x = Oxygen, y = fraction, colour = Oxygen)
) +
  geom_point(
    position = position_jitter(width = 0.09, height = 0),
    size = 2.0,
    alpha = 0.78
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.48,
    linewidth = 0.65,
    colour = ink
  ) +
  facet_grid(metric_label ~ Age) +
  scale_colour_manual(values = c(
    "Normoxia" = air_colour,
    "Hyperoxia" = hyperoxia_colour
  )) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(
    title = "B. Animal-level capillary composition",
    subtitle = "Points are individual animals; horizontal bars are group means",
    x = NULL,
    y = "Fraction of endothelial cells",
    colour = NULL
  ) +
  theme_research +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "none"
  )

test_metric_order <- c(
  "endothelial_fraction_all_cells",
  "Cap_fraction_of_endothelium",
  "Cap-a_fraction_of_endothelium",
  "Art_fraction_of_endothelium",
  "Vein_fraction_of_endothelium",
  "Lymph_fraction_of_endothelium"
)
exact_composition_labels <- c(
  "endothelial_fraction_all_cells" = "Endothelial / all cells",
  "Cap_fraction_of_endothelium" = "gCap-like",
  "Cap-a_fraction_of_endothelium" = "aCap-like",
  "Art_fraction_of_endothelium" = "Arterial",
  "Vein_fraction_of_endothelium" = "Venous",
  "Lymph_fraction_of_endothelium" = "Lymphatic"
)
composition_tests[, metric_label := factor(
  exact_composition_labels[metric],
  levels = rev(exact_composition_labels[test_metric_order])
)]
composition_tests[, Age := factor(Age, levels = c("P3", "P7", "P14"))]
composition_tests[, fdr_status := fifelse(
  BH_q_all_18_logit_tests < 0.05,
  "FDR < 0.05",
  "FDR >= 0.05"
)]
composition_tests[, `:=`(
  effect_pp = 100 * raw_fraction_difference_hyperoxia_minus_normoxia,
  ci_low_pp = 100 * raw_difference_ci95_low,
  ci_high_pp = 100 * raw_difference_ci95_high
)]

p3c <- ggplot(
  composition_tests,
  aes(x = effect_pp, y = metric_label, colour = fdr_status)
) +
  geom_vline(xintercept = 0, colour = muted, linewidth = 0.5) +
  geom_errorbar(
    aes(xmin = ci_low_pp, xmax = ci_high_pp),
    orientation = "y",
    width = 0.17,
    linewidth = 0.75
  ) +
  geom_point(size = 2.7) +
  facet_grid(. ~ Age) +
  scale_colour_manual(values = c(
    "FDR < 0.05" = hyperoxia_colour,
    "FDR >= 0.05" = "#8B969E"
  )) +
  labs(
    title = "C. Hyperoxia-associated composition differences",
    subtitle = "Raw percentage-point difference; 10,000-repeat within-group bootstrap 95% interval; colour uses exact-permutation BH q across 18 tests",
    x = "Hyperoxia minus normoxia (percentage points)",
    y = NULL,
    colour = NULL
  ) +
  theme_research +
  theme(legend.position = "top")

f3_top <- plot_grid(p3a, p3b, ncol = 2, rel_widths = c(0.92, 1.08))
f3_body <- plot_grid(f3_top, p3c, ncol = 1, rel_heights = c(1.05, 0.95))
f3_final <- add_figure_header(
  f3_body,
  "Figure 3. Hyperoxia is associated with shifts in recovered endothelial composition",
  "Animal-level GSE151974 analyses show the largest compositional association at postnatal day 14; recovery and state abundance are not interpreted causally."
)
save_main_figure(f3_final, "Figure3_GSE151974_composition", 16, 11.8)
fwrite(composition_means, file.path(plot_data_dir, "Figure3_composition_means.tsv"), sep = "\t")
fwrite(composition_long, file.path(plot_data_dir, "Figure3_animal_composition.tsv"), sep = "\t")
fwrite(composition_tests, file.path(plot_data_dir, "Figure3_composition_tests.tsv"), sep = "\t")

# -------------------------------------------------------------------------
# Figure 4. Cross-dataset gene and pathway replication
# -------------------------------------------------------------------------

scatter_effects <- rbindlist(lapply(c("Cap", "Cap-a"), function(cell_type) {
  pb <- pseudobulk[
    model == "age_by_oxygen_P14_contrast" & CellType == cell_type,
    .(
      gene_upper = gene_canonical,
      pseudobulk_logFC = logFC,
      pseudobulk_FDR = FDR
    )
  ]
  z <- merge(
    bulk_de[, .(
      gene_upper = gene_canonical,
      bulk_log2FC = log2FoldChange_mle,
      bulk_FDR = padj
    )],
    pb,
    by = "gene_upper"
  )
  z[, CellType := cell_type]
  z
}))
scatter_effects[, both_fdr05 := bulk_FDR < 0.05 & pseudobulk_FDR < 0.05]
scatter_effects[, core_gene := FALSE]

replication_annotations <- replication_summary[
  CellType %in% c("Cap", "Cap-a"),
  .(
    CellType,
    label = sprintf(
      "rho=%.2f | %d/%d overlap genes concordant",
      spearman_all_shared_gene_effects,
      direction_concordant_overlap_genes,
      significant_overlap_genes
    )
  )
]

p4a <- ggplot(
  scatter_effects,
  aes(x = bulk_log2FC, y = pseudobulk_logFC)
) +
  geom_hline(yintercept = 0, colour = grid_colour, linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = grid_colour, linewidth = 0.35) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#AAB2B8", linewidth = 0.5) +
  geom_point(colour = "#7F8B94", size = 0.55, alpha = 0.17) +
  geom_point(
    data = scatter_effects[both_fdr05 == TRUE],
    aes(colour = pseudobulk_logFC > 0),
    size = 0.9,
    alpha = 0.52
  ) +
  geom_text(
    data = replication_annotations,
    aes(x = -2.8, y = 4.65, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    colour = muted,
    size = 2.65
  ) +
  facet_grid(. ~ CellType) +
  scale_colour_manual(values = c("TRUE" = up_colour, "FALSE" = down_colour), guide = "none") +
  coord_cartesian(xlim = c(-3.1, 6.0), ylim = c(-3.7, 5.1)) +
  labs(
    title = "A. Cross-platform gene-effect agreement",
    subtitle = "Bulk DESeq2 versus animal x subtype pseudobulk edgeR (categorical age x oxygen P14 contrast)",
    x = "GSE216046 bulk log2 fold-change",
    y = "GSE151974 pseudobulk log2 fold-change"
  ) +
  theme_research +
  theme(legend.position = "none")

replication_summary[, CellType := factor(
  CellType,
  levels = c("Lymph", "Art", "Vein", "Cap-a", "Cap")
)]
replication_summary[, concordance_label := sprintf(
  "%d/%d",
  direction_concordant_overlap_genes,
  significant_overlap_genes
)]

p4b <- ggplot(
  replication_summary,
  aes(
    x = spearman_all_shared_gene_effects,
    y = CellType,
    size = significant_overlap_genes,
    fill = direction_concordance_fraction
  )
) +
  geom_segment(
    aes(x = 0, xend = spearman_all_shared_gene_effects, yend = CellType),
    colour = "#BEC6CC",
    linewidth = 0.8
  ) +
  geom_point(shape = 21, colour = ink, stroke = 0.7) +
  geom_text(
    aes(label = concordance_label),
    nudge_x = 0.055,
    hjust = 0,
    colour = ink,
    size = 2.6
  ) +
  scale_fill_gradient(low = "#DDE7EE", high = green_colour, limits = c(0.8, 1)) +
  scale_size_continuous(range = c(4, 10)) +
  scale_x_continuous(limits = c(0, 0.83), expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "B. Replication strength by subtype",
    subtitle = "Labels: concordant / overlap genes; darker fill, higher concordance",
    x = "Spearman correlation across all shared genes",
    y = NULL,
    size = "Overlap genes",
    fill = "Direction\nconcordance"
  ) +
  theme_research +
  theme(legend.position = "none")

selected_pathways <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
  "REACTOME_COMPLEX_I_BIOGENESIS",
  "REACTOME_NUCLEAR_EVENTS_MEDIATED_BY_NFE2L2",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX"
)

pathway_heatmap <- rbind(
  fgsea_all[
    dataset == "GSE216046" & pathway %in% selected_pathways,
    .(pathway, evidence_column = "Bulk discovery", NES, padj)
  ],
  fgsea_all[
    dataset == "GSE151974"
      & model == "age_by_oxygen_P14_contrast"
      & pathway %in% selected_pathways,
    .(pathway, evidence_column = cell_type, NES, padj)
  ]
)
pathway_heatmap[, pathway_label := clean_pathway_label(pathway)]
pathway_heatmap[, pathway_label := factor(
  pathway_label,
  levels = rev(clean_pathway_label(selected_pathways))
)]
pathway_heatmap[, evidence_column := factor(
  evidence_column,
  levels = c("Bulk discovery", "Cap", "Cap-a", "Art", "Vein", "Lymph")
)]
pathway_heatmap[, significance := fifelse(padj < 0.05, "*", "")]
pathway_heatmap[, nes_fill := pmax(-3, pmin(3, NES))]

p4c <- ggplot(
  pathway_heatmap,
  aes(x = evidence_column, y = pathway_label, fill = nes_fill)
) +
  geom_tile(colour = "white", linewidth = 0.65) +
  geom_text(
    aes(label = sprintf("%.2f%s", NES, significance)),
    size = 2.55,
    colour = ink
  ) +
  scale_fill_gradient2(
    low = down_colour,
    mid = "white",
    high = up_colour,
    midpoint = 0,
    limits = c(-3, 3)
  ) +
  labs(
    title = "C. Direction-aware pathway replication",
    subtitle = "NES values; *FDR<0.05 within each dataset, subtype and collection",
    x = NULL,
    y = NULL,
    fill = "NES"
  ) +
  theme_research +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "right"
  )

f4_top <- plot_grid(p4a, p4b, ncol = 2, rel_widths = c(1.38, 0.62))
f4_body <- plot_grid(f4_top, p4c, ncol = 1, rel_heights = c(0.92, 1.08))
f4_final <- add_figure_header(
  f4_body,
  "Figure 4. Capillary endothelium shows the strongest gene and pathway replication",
  "Replication requires independent mouse-dataset significance, |observed log2FC|>=1, and direction concordance; the P14 interaction contrast is primary."
)
save_main_figure(f4_final, "Figure4_mouse_replication", 16.5, 12.5)
fwrite(scatter_effects, file.path(plot_data_dir, "Figure4_gene_effect_scatter.tsv"), sep = "\t")
fwrite(replication_summary, file.path(plot_data_dir, "Figure4_replication_summary.tsv"), sep = "\t")
fwrite(pathway_heatmap, file.path(plot_data_dir, "Figure4_pathway_heatmap.tsv"), sep = "\t")

# -------------------------------------------------------------------------
# Figure 5. Human projection (already independently generated and audited)
# -------------------------------------------------------------------------

human_source_png <- file.path(
  project_root,
  "08_figures",
  "GSE275938_human_endothelial_subtype_sensitivity.png"
)
human_source_pdf <- file.path(
  project_root,
  "08_figures",
  "GSE275938_human_endothelial_subtype_sensitivity.pdf"
)
stopifnot(file.exists(human_source_png), file.exists(human_source_pdf))
file.copy(
  human_source_png,
  file.path(figure_dir, "Figure5_human_subtype_projection.png"),
  overwrite = TRUE
)
file.copy(
  human_source_pdf,
  file.path(figure_dir, "Figure5_human_subtype_projection.pdf"),
  overwrite = TRUE
)
fwrite(
  human_compact,
  file.path(plot_data_dir, "Figure5_human_subtype_compact_evidence.tsv"),
  sep = "\t"
)

# -------------------------------------------------------------------------
# Figure 6. Exploratory post-projection subset and evidence annotations
# -------------------------------------------------------------------------

core_matrix <- fread(file.path(
  result_root,
  "human_endothelial_subtypes",
  "GSE275938_exploratory_subset_donor_counts_and_log2CPM.tsv"
))[
  subtype %in% c("gCap", "aCap")
]
core_matrix[, gene_upper := toupper(gene)]
core_matrix[, expression_mean := mean(log2CPM_TMM), by = .(gene_upper, subtype)]
core_matrix[, expression_sd := sd(log2CPM_TMM), by = .(gene_upper, subtype)]
core_matrix[, expression_z := fifelse(
  expression_sd > 0,
  (log2CPM_TMM - expression_mean) / expression_sd,
  0
)]
core_matrix[, evidence_column := paste0(
  subtype,
  " | ",
  fifelse(
    sample == "BPD+PH 1",
    "PH1",
    fifelse(
      sample == "BPD+PH 2",
      "PH2",
      gsub(" ", "", sample)
    )
  )
)]
evidence_columns <- c(
  "gCap | BPD1",
  "gCap | BPD2",
  "gCap | PH1",
  "gCap | PH2",
  "aCap | BPD1",
  "aCap | BPD2",
  "aCap | PH1",
  "aCap | PH2"
)
core_matrix[, gene_upper := factor(gene_upper, levels = rev(core_labels))]
core_matrix[, evidence_column := factor(
  evidence_column,
  levels = evidence_columns
)]
core_matrix[, expression_label := sprintf("%.1f", log2CPM_TMM)]

p6a <- ggplot(
  core_matrix,
  aes(x = evidence_column, y = gene_upper, fill = expression_z)
) +
  geom_tile(colour = "white", linewidth = 0.75) +
  geom_text(aes(label = expression_label), colour = ink, size = 2.25) +
  scale_fill_gradient2(
    low = down_colour,
    mid = "white",
    high = up_colour,
    midpoint = 0,
    limits = c(-1.5, 1.5),
    oob = scales::squish
  ) +
  labs(
    title = sprintf(
      "A. Donor-level expression of %d post-projection genes",
      length(core_labels)
    ),
    subtitle = "Numbers are TMM log2CPM; colour is standardized within each gene and subtype. Groups: BPD and BPD+PH.",
    x = NULL,
    y = NULL,
    fill = "Within-gene\nexpression z"
  ) +
  theme_research +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 27, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  )

mechanism_boxes <- data.table(
  id = c("mouse", "respiration", "stress", "ecm", "composition", "human"),
  xmin = c(0.02, 0.35, 0.68, 0.02, 0.35, 0.68),
  xmax = c(0.32, 0.65, 0.98, 0.32, 0.65, 0.98),
  ymin = c(0.56, 0.56, 0.56, 0.08, 0.08, 0.08),
  ymax = c(0.94, 0.94, 0.94, 0.46, 0.46, 0.46),
  fill = c("#E4EEF5", "#E4EEF5", "#F6E7E3", "#EEE7F4", "#E8F0EA", "#FFF2EB"),
  border = c(air_colour, air_colour, up_colour, purple_colour, green_colour, hyperoxia_colour),
  label = c(
    "Mouse replication layer\nP14 age x oxygen contrast\n198 directional genes\n33 genes in >=3 subtypes",
    "Respiration-annotated sets\nOxidative phosphorylation\nand electron transport\nDirectional replication",
    "Stress-program annotations\np53, ROS and apoptosis\nTranscript-level evidence\nNo causal mediation claim",
    "Matrix-program annotations\nECM organization and turnover\nTranscript-level evidence\nNo structural validation",
    "Recovered composition\nP14 CLR permutation P=0.0022\nAssociation with recovered states\nNot causal remodeling",
    "Human projection\nBPD n=2 | BPD+PH n=2\n12.8-fold recovery imbalance\nExploratory support only"
  )
)

p6b <- ggplot() +
  geom_rect(
    data = mechanism_boxes,
    aes(
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
      fill = fill, colour = border
    ),
    linewidth = 1.0,
    show.legend = FALSE
  ) +
  scale_fill_identity() +
  scale_colour_identity() +
  geom_text(
    data = mechanism_boxes,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    colour = ink,
    size = 3.05,
    lineheight = 1.03,
    fontface = ifelse(
      mechanism_boxes$id %in% c("mouse", "composition", "human"),
      "bold",
      "plain"
    )
  ) +
  annotate(
    "text",
    x = 0.5,
    y = 0.005,
    label = "Independent evidence layers are juxtaposed without causal edges; pathway terms are transcription-program annotations.",
    size = 3.0,
    colour = muted,
    fontface = "italic"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(
    title = "B. Evidence-layer annotations without causal links",
    subtitle = "The P14 interaction model is primary; all-age, cell-threshold and effect-threshold analyses are sensitivities."
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", colour = ink, size = 11, margin = margin(b = 3)),
    plot.subtitle = element_text(colour = muted, size = 8.8, margin = margin(b = 5)),
    plot.margin = margin(8, 10, 8, 8)
  )

f6_body <- plot_grid(p6a, p6b, ncol = 1, rel_heights = c(1.05, 0.95))
f6_final <- add_figure_header(
  f6_body,
  "Figure 6. Exploratory human-concordant subset and evidence annotations",
  "The post-projection subset is hypothesis-generating; donor-level human expression is shown explicitly and mechanistic language is avoided."
)
save_main_figure(f6_final, "Figure6_core_and_model", 16.5, 12)
fwrite(core_matrix, file.path(plot_data_dir, "Figure6_core_gene_matrix.tsv"), sep = "\t")
fwrite(mechanism_boxes, file.path(plot_data_dir, "Figure6_mechanism_boxes.tsv"), sep = "\t")

writeLines(
  capture.output(sessionInfo()),
  file.path(figure_dir, "manuscript_main_figures_sessionInfo.txt")
)

cat(
  sprintf(
    "Created manuscript figures in %s\nCreated plotting data in %s\n",
    figure_dir,
    plot_data_dir
  )
)
