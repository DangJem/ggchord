# ============================================================================
# make_readme_plots.R
# Regenerates every image referenced by README.md / README-Hans.md into
# man/figures/ (the location required by pkgdown for README images).
#
# The README examples use the built-in datasets (seq_data_example,
# ribbon_data_example, gene_data_example), so the code blocks in the README
# can be copy-pasted and run directly.
#
# Usage:
#   Rscript examples/make_readme_plots.R
# ============================================================================

suppressPackageStartupMessages(pkgload::load_all("."))

data(seq_data_example)
data(ribbon_data_example)
data(gene_data_example)

outdir <- file.path("man", "figures")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, width = 800, height = 800, res = 100, bg = NULL) {
  if (is.null(bg)) {
    grDevices::png(file.path(outdir, name), width = width, height = height,
                   res = res)
  } else {
    grDevices::png(file.path(outdir, name), width = width, height = height,
                   res = res, bg = bg)
  }
  suppressMessages(suppressWarnings(print(p)))
  grDevices::dev.off()
  cat("saved:", name, "\n")
}

# ---------------------------------------------------------------------------
# 1. Quick start: all defaults
# ---------------------------------------------------------------------------
save_plot(
  ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() + geom_ribbon() + geom_gene() + geom_axis(),
  "combined_default.png")

# ---------------------------------------------------------------------------
# 2. Sequences only
# ---------------------------------------------------------------------------
save_plot(
  ggchord(seq_data_example) + geom_seq(),
  "seq_only_default.png")

save_plot(
  ggchord(seq_data_example) +
    geom_seq(
      seq_order = c("MT118296.1", "OR222515.1", "MT108731.1", "OQ646790.1"),
      seq_orientation = c(1, -1, 1, -1),
      seq_curvature = c(0, 2, -2, 6),
      seq_colors = c("steelblue", "orange", "pink", "yellow")
    ),
  "seq_only_custom.png")

# ---------------------------------------------------------------------------
# 3. Alignment ribbons: color schemes and outlines
# ---------------------------------------------------------------------------
save_plot(
  ggchord(seq_data_example, ribbon_data_example) +
    geom_seq() + geom_ribbon(),
  "ribbon_pident.png")

save_plot(
  ggchord(seq_data_example, ribbon_data_example) +
    geom_seq() + geom_ribbon(ribbon_color_scheme = "query"),
  "ribbon_query.png")

save_plot(
  ggchord(seq_data_example, ribbon_data_example) +
    geom_seq() + geom_ribbon(ribbon_color_scheme = "subject"),
  "ribbon_subject.png")

save_plot(
  ggchord(seq_data_example, ribbon_data_example) +
    geom_seq() +
    geom_ribbon(ribbon_color_scheme = "single", ribbon_colors = "orange"),
  "ribbon_single.png")

# ---------------------------------------------------------------------------
# 4. Gene annotations
# ---------------------------------------------------------------------------
save_plot(
  ggchord(seq_data_example, gene_data = gene_data_example) +
    geom_seq() +
    geom_gene() +
    geom_gene_label_repel(gene_label_wrap = 15, max_overlaps = 5),
  "gene_repel.png")
save_plot(
  ggchord(seq_data_example, gene_data = gene_data_example) +
    geom_seq() +
    geom_gene() +
    geom_gene_label_repel(
      gene_label_orientation = "horizontal",
      gene_label_segment = "elbow",
      gene_label_wrap = 15
    ),
  "gene_repel_elbow.png")
save_plot(
  ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() +
    geom_ribbon() +
    geom_gene() +
    geom_gene_label_repel(
      gene_label_orientation = "horizontal",
      gene_label_segment = "elbow",
      gene_label_side = "outside",
      gene_label_wrap = 15
    ),
  "gene_repel_outside.png")
save_plot(
  ggchord(seq_data_example, gene_data = gene_data_example) +
    geom_seq() + geom_gene(),
  "gene_strand.png")

save_plot(
  ggchord(seq_data_example, gene_data = gene_data_example) +
    geom_seq() +
    geom_gene(gene_color_scheme = "manual") +
    geom_gene_label(
      gene_label_rotation = 45,
      gene_label_radial_offset = 0.1
    ),
  "gene_manual_label.png")

# ---------------------------------------------------------------------------
# 5. Axes and sequence labels
# ---------------------------------------------------------------------------
save_plot(
  ggchord(seq_data_example) +
    geom_seq() +
    geom_axis(axis_tick_major_length = 0.03,
              axis_tick_minor_length = 0.015,
              axis_label_size = 2.5) +
    geom_seq_label(seq_label_radius = 1.2),
  "axis_seq_label.png")
save_plot(
  ggchord(seq_data_example, rotation = 30) +
    geom_seq() +
    geom_seq_label(seq_label_radius = 1.15,
                   seq_label_orientation = "horizontal",
                   seq_label_size = 3.5,
                   colour = "#2563EB"),
  "seq_label_horizontal.png")

# ---------------------------------------------------------------------------
# 6. Two-sequence and three-sequence comparisons (subset of the example data)
# ---------------------------------------------------------------------------
seq2 <- seq_data_example[seq_data_example$seq_id %in% c("MT108731.1", "MT118296.1"), ]
ribbon2 <- ribbon_data_example[
  ribbon_data_example$qaccver %in% seq2$seq_id &
  ribbon_data_example$saccver %in% seq2$seq_id, ]
gene2 <- gene_data_example[gene_data_example$seq_id %in% seq2$seq_id, ]
save_plot(
  ggchord(seq2, ribbon2, gene2) + geom_seq() + geom_ribbon() + geom_gene() + geom_axis(),
  "example_seq2.png")

seq3 <- seq_data_example[seq_data_example$seq_id %in% c("MT108731.1", "MT118296.1", "OQ646790.1"), ]
ribbon3 <- ribbon_data_example[
  ribbon_data_example$qaccver %in% seq3$seq_id &
  ribbon_data_example$saccver %in% seq3$seq_id, ]
gene3 <- gene_data_example[gene_data_example$seq_id %in% seq3$seq_id, ]
save_plot(
  ggchord(seq3, ribbon3, gene3) + geom_seq() + geom_ribbon() + geom_gene() + geom_axis(),
  "example_seq3.png")

# ---------------------------------------------------------------------------
# 7. Full customization
# ---------------------------------------------------------------------------
# A cover-ready showcase: every layer contributes its distinctive feature --
# per-sequence radii/orientation/curvature, pident ribbons, strand-colored
# genes, repelled labels (horizontal + elbow + moved outside the arcs),
# sequence labels, axes, a custom theme, and a curated palette.
save_plot(
  ggchord(
    seq_data = seq_data_example,
    ribbon_data = ribbon_data_example,
    gene_data = gene_data_example,
    title = "ggchord",
    rotation = 45,
    panel_margin = list(t = 1.5, r = 0.6, b = 0.6, l = 0.6)
  ) +
    labs(subtitle = "Layered multi-sequence chord diagrams for ggplot2") +
    geom_seq(
      seq_radius = c(3.3, 2.5, 1.8, 1.25),
      seq_orientation = c(-1, -1, 1, -1),
      seq_curvature = c(0.8, 1, 1.4, 1),
      seq_gap = 0.035,
      seq_colors = c(
        "MT108731.1" = "#E76F51",
        "MT118296.1" = "#264653",
        "OQ646790.1" = "#2A9D8F",
        "OR222515.1" = "#D9A62E"
      ),
      linewidth = 1.6
    ) +
    geom_ribbon(
      ribbon_color_scheme = "pident",
      ribbon_gap = 0.12,
      ribbon_alpha = 0.45,
      ribbon_outline_color = "#FBF9F6",
      ribbon_outline_width = 0.03
    ) +
    geom_gene(
      gene_offset = 0.1,
      gene_width = 0.06,
      gene_colors = c("+" = "#4C6EF5", "-" = "#F06595")
    ) +
    geom_gene_label_repel(
      gene_label_orientation = "horizontal",
      gene_label_segment = "elbow",
      gene_label_side = "outside",
      gene_label_wrap = 12,
      gene_label_size = 2,
      seed = 42
    ) +
    geom_seq_label(
      seq_label_radius = 1.1,
      seq_label_size = 3.4,
      colour = "#52525B"
    ) +
    geom_axis(
      axis_gap = 0.07,
      axis_tick_major_number = 4,
      axis_tick_major_length = 0.025,
      axis_tick_minor_number = 4,
      axis_tick_minor_length = 0.012,
      axis_label_size = 1.8
    ) +
    theme(
      plot.background = element_rect(fill = "#FBF9F6", colour = NA),
      panel.background = element_rect(fill = "#FBF9F6", colour = NA),
      plot.title = element_text(
        size = 26, face = "bold", colour = "#1F2937",
        hjust = 0.5, margin = margin(t = 10, b = 2)
      ),
      plot.subtitle = element_text(
        size = 12, colour = "#6B7280",
        hjust = 0.5, margin = margin(b = 12)
      ),
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold", colour = "#374151"),
      legend.text = element_text(size = 8, colour = "#4B5563"),
      legend.key.size = unit(0.7, "cm")
    ),
  "combined_fine.png",
  width = 1600, height = 1600, res = 200, bg = "#FBF9F6")

# ---------------------------------------------------------------------------
# 7. Themes and scales added with "+"
# ---------------------------------------------------------------------------
save_plot(
  ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() + geom_ribbon() + geom_gene() + geom_axis() +
    theme(panel.background = element_rect(fill = "grey95")) +
    scale_color_manual(
      values = c("MT108731.1" = "#E41A1C",
                 "MT118296.1" = "#377EB8",
                 "OQ646790.1" = "#4DAF4A",
                 "OR222515.1" = "#984EA3")
    ),
  "theme_custom.png")

# ---------------------------------------------------------------------------
# 8b. Legends grouped with theme()
# ---------------------------------------------------------------------------
save_plot(
  ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq(legend_position = NULL) +
    geom_ribbon(legend_position = NULL) +
    geom_gene(legend_position = NULL) +
    geom_axis() +
    theme(legend.position = "bottom", legend.box = "horizontal"),
  "legend_bottom.png")

cat("ALL PLOTS DONE\n")
