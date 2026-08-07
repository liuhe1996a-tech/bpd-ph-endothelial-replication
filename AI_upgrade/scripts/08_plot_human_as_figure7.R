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
})

result_dir <- file.path(
  project_root, "07_results", "human_endothelial_subtypes"
)
figure_dir <- file.path(project_root, "08_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

ink <- "#27313A"
blue <- "#2C5F8A"
orange <- "#D97932"
light_blue <- "#DCE8F1"
light_orange <- "#F4E2D3"
grid_colour <- "#D9DEE3"

sample_order <- c("BPD 1", "BPD 2", "BPD+PH 1", "BPD+PH 2")
subtype_order <- c(
  "All Endothelial",
  "gCap",
  "aCap",
  "Arterial EC",
  "Lymphatic",
  "Pulmonary venous EC"
)
subtype_labels <- c(
  "All Endothelial" = "All EC",
  "gCap" = "gCap",
  "aCap" = "aCap",
  "Arterial EC" = "Arterial EC",
  "Lymphatic" = "Lymphatic EC",
  "Pulmonary venous EC" = "Pulmonary venous EC",
  "Systemic venous EC" = "Systemic venous EC",
  "abCap" = "abCap"
)

theme_research <- theme_minimal(base_size = 10) +
  theme(
    text = element_text(colour = ink),
    axis.text = element_text(colour = ink),
    axis.title = element_text(colour = ink),
    plot.title = element_text(face = "bold", colour = ink, size = 11),
    plot.subtitle = element_text(colour = "#5A6670", size = 9),
    panel.grid.major = element_line(colour = grid_colour, linewidth = 0.35),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", colour = ink),
    plot.margin = margin(8, 12, 8, 8)
  )

coverage <- fread(
  file.path(
    result_dir,
    "GSE275938_disease_endothelial_subtype_cell_coverage.tsv"
  )
)[subtype != "All Endothelial"]
coverage_long <- melt(
  coverage,
  id.vars = "subtype",
  measure.vars = sample_order,
  variable.name = "sample",
  value.name = "cells"
)
coverage_long[, sample := factor(sample, levels = sample_order)]
coverage_long[, subtype := factor(
  subtype,
  levels = rev(c(
    "gCap", "aCap", "abCap", "Arterial EC",
    "Pulmonary venous EC", "Systemic venous EC", "Lymphatic"
  )),
  labels = rev(subtype_labels[c(
    "gCap", "aCap", "abCap", "Arterial EC",
    "Pulmonary venous EC", "Systemic venous EC", "Lymphatic"
  )])
)]

panel_a <- ggplot(
  coverage_long,
  aes(x = sample, y = subtype, fill = log10(cells + 1))
) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = cells), colour = ink, size = 3.0) +
  scale_fill_gradient(
    low = light_blue,
    high = blue,
    name = expression(log[10](cells + 1))
  ) +
  labs(
    title = "A. Endothelial-cell recovery by donor and subtype",
    subtitle = "Exact cell counts; BPD and BPD+PH each contain two donors",
    x = NULL,
    y = NULL
  ) +
  theme_research +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "right"
  )

donor_scores <- fread(
  file.path(
    result_dir,
    "GSE275938_subtype_locked_signature_donor_scores.tsv"
  )
)[
  normalization == "TMM"
    & signature == "mouse_subtypes_ge3"
    & subtype %in% subtype_order
]
eligible_labels <- unique(donor_scores[, .(
  subtype,
  genes_eligible,
  genes_requested
)])
eligible_labels[, display := sprintf(
  "%s (%d/%d genes)",
  subtype_labels[subtype],
  genes_eligible,
  genes_requested
)]
display_levels <- eligible_labels[
  match(subtype_order, subtype),
  display
]
donor_scores <- merge(
  donor_scores,
  eligible_labels[, .(subtype, display)],
  by = "subtype",
  all.x = TRUE
)
donor_scores[, display := factor(
  display,
  levels = rev(display_levels)
)]
donor_scores[, condition := factor(
  condition,
  levels = c("BPD", "BPD+PH")
)]

panel_b <- ggplot(
  donor_scores,
  aes(
    x = oriented_signature_score,
    y = display,
    colour = condition,
    shape = condition
  )
) +
  geom_vline(xintercept = 0, colour = "#808991", linewidth = 0.5) +
  geom_point(
    size = 3.1,
    stroke = 0.9,
    position = position_dodge(width = 0.42)
  ) +
  scale_colour_manual(values = c("BPD" = blue, "BPD+PH" = orange)) +
  scale_shape_manual(values = c("BPD" = 16, "BPD+PH" = 17)) +
  labs(
    title = "B. Donor-level mouse-signature scores",
    subtitle = "TMM scores; labels give eligible / 33 one-to-one multi-subtype genes",
    x = "Oriented signature score",
    y = NULL
  ) +
  theme_research

downsampling <- fread(
  file.path(
    result_dir,
    "GSE275938_equal_cell_downsampling_summary.tsv"
  )
)[
  signature == "mouse_subtypes_ge3"
    & subtype %in% subtype_order
]
downsampling[, method := "Equal cell count"]
downsampling[, target_label := paste0(cells_per_donor, " cells")]
rarefaction <- fread(file.path(
  result_dir,
  "GSE275938_library_depth_rarefaction_summary.tsv"
))[
  signature == "mouse_subtypes_ge3"
    & subtype %in% subtype_order
]
rarefaction[, method := "Equal total UMI depth"]
rarefaction[, target_label := paste0(
  format(target_total_umi_per_donor, big.mark = ","),
  " UMI"
)]
sampling_sensitivity <- rbindlist(
  list(
    downsampling[, .(
      subtype,
      method,
      target_label,
      median_effect,
      sampling_p2_5,
      sampling_p97_5
    )],
    rarefaction[, .(
      subtype,
      method,
      target_label,
      median_effect,
      sampling_p2_5,
      sampling_p97_5
    )]
  )
)
sampling_sensitivity[, subtype := factor(
  subtype,
  levels = rev(subtype_order),
  labels = rev(subtype_labels[subtype_order])
)]

panel_c <- ggplot(
  sampling_sensitivity,
  aes(x = median_effect, y = subtype, colour = method, fill = method)
) +
  geom_vline(xintercept = 0, colour = "#808991", linewidth = 0.5) +
  geom_errorbar(
    aes(xmin = sampling_p2_5, xmax = sampling_p97_5),
    orientation = "y",
    width = 0.18,
    linewidth = 0.9,
    position = position_dodge(width = 0.45)
  ) +
  geom_point(
    shape = 21,
    size = 3.2,
    stroke = 0.9,
    position = position_dodge(width = 0.45)
  ) +
  scale_colour_manual(values = c(
    "Equal cell count" = orange,
    "Equal total UMI depth" = blue
  )) +
  scale_fill_manual(values = c(
    "Equal cell count" = light_orange,
    "Equal total UMI depth" = light_blue
  )) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  labs(
    title = "C. Cell-count and UMI-depth sensitivity",
    subtitle = "Median and 2.5th–97.5th sampling percentiles across 500 repeats per method",
    x = "BPD+PH minus BPD signature-score effect",
    y = NULL
  ) +
  theme_research

random_null <- fread(
  file.path(
    project_root,
    "07_results",
    "sensitivity_analyses",
    "GSE275938_six_population_random_set_results.tsv"
  )
)[
  signature == "mouse_subtypes_ge3"
    & subtype %in% subtype_order
]
random_null[, subtype := factor(
  subtype,
  levels = rev(subtype_order),
  labels = rev(subtype_labels[subtype_order])
)]
random_null[, p_label := sprintf(
  "BH q=%.3f",
  empirical_q_greater_BH_within_signature
)]

panel_d <- ggplot(
  random_null,
  aes(y = subtype)
) +
  geom_vline(xintercept = 0, colour = "#808991", linewidth = 0.5) +
  geom_errorbar(
    aes(xmin = null_p2_5, xmax = null_p97_5),
    orientation = "y",
    width = 0.18,
    colour = "#88939C",
    linewidth = 1.0
  ) +
  geom_point(
    aes(x = null_median),
    shape = 21,
    size = 2.5,
    stroke = 0.7,
    colour = "#69747D",
    fill = "white"
  ) +
  geom_point(
    aes(x = observed_effect),
    shape = 18,
    size = 3.4,
    colour = blue
  ) +
  geom_text(
    aes(x = pmax(observed_effect, null_p97_5), label = p_label),
    hjust = -0.08,
    colour = "#5A6670",
    size = 2.6
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.35))) +
  labs(
    title = "D. Orthology- and expression-matched specificity",
    subtitle = "Mouse signature versus 2,000 matched random sets; BH over six populations",
    x = "BPD+PH minus BPD oriented effect",
    y = NULL
  ) +
  theme_research

combined <- plot_grid(
  panel_a,
  panel_b,
  panel_c,
  panel_d,
  ncol = 2,
  rel_widths = c(1.03, 1),
  rel_heights = c(1, 1),
  align = "hv"
)
title <- ggdraw() +
  draw_label(
    "Figure 7. Exploratory human GSE275938 subtype projection",
    x = 0.01,
    hjust = 0,
    fontface = "bold",
    size = 15,
    colour = ink
  ) +
  draw_label(
    "Donor is the unit (n=2+2); diagnosis is confounded with cell recovery and sequencing depth",
    x = 0.01,
    y = 0.18,
    hjust = 0,
    size = 9.5,
    colour = "#5A6670"
  )
final_plot <- plot_grid(
  title,
  combined,
  ncol = 1,
  rel_heights = c(0.08, 1)
)

ggsave(
  file.path(
    figure_dir,
    "GSE275938_human_endothelial_subtype_sensitivity.png"
  ),
  final_plot,
  width = 14,
  height = 10,
  dpi = 320,
  bg = "white"
)
ggsave(
  file.path(
    figure_dir,
    "GSE275938_human_endothelial_subtype_sensitivity.pdf"
  ),
  final_plot,
  width = 14,
  height = 10,
  bg = "white"
)
writeLines(
  capture.output(sessionInfo()),
  file.path(
    result_dir,
    "GSE275938_human_subtype_figure_sessionInfo.txt"
  )
)
