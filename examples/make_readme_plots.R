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

save_plot <- function(p, name, width = 800, height = 800, res = 100) {
  grDevices::png(file.path(outdir, name), width = width, height = height, res = res)
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
    geom_seq() + geom_gene(),
  "gene_strand.png")

save_plot(
  ggchord(seq_data_example, gene_data = gene_data_example) +
    geom_seq() +
    geom_gene(
      gene_color_scheme = "manual",
      gene_label_show = TRUE,
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
save_plot(
  ggchord(
    seq_data = seq_data_example,
    ribbon_data = ribbon_data_example,
    gene_data = gene_data_example,
    title = "Multi-sequence Chord Diagram with Gene Annotations",
    rotation = 45
  ) +
    geom_seq(
      seq_radius = c(3, 2, 2, 1),
      seq_orientation = c(-1, -1, -1, 1),
      seq_curvature = c(0, 1, -1, 1.5),
      seq_gap = 0.03
    ) +
    geom_ribbon(ribbon_color_scheme = "pident", ribbon_gap = 0.1) +
    geom_gene(
      gene_offset = list(
        c("+" = 0.2, "-" = -0.2),
        c("+" = 0.2, "-" = -0.2),
        c("+" = 0.2, "-" = 0),
        c("+" = 0.2, "-" = 0.1)
      ),
      gene_width = 0.08,
      gene_label_show = TRUE,
      gene_label_rotation = list(
        c("+" = 45, "-" = -45),
        c("+" = 30, "-" = -30),
        c("+" = 15, "-" = -15),
        c("+" = 0, "-" = 0)
      )
    ) +
    geom_axis(
      axis_gap = 0,
      axis_tick_major_length = 0.03,
      axis_label_size = 2
    ),
  "combined_fine.png")

# ---------------------------------------------------------------------------
# 7. Themes and scales added with "+"
# ---------------------------------------------------------------------------
save_plot(
  ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() + geom_ribbon() + geom_gene() + geom_axis() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.background = element_rect(fill = "grey95")
    ) +
    scale_color_manual(
      values = c("MT108731.1" = "#E41A1C",
                 "MT118296.1" = "#377EB8",
                 "OQ646790.1" = "#4DAF4A",
                 "OR222515.1" = "#984EA3")
    ),
  "theme_custom.png")

cat("ALL PLOTS DONE\n")
