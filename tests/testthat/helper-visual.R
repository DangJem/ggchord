# Helper functions for the lightweight visual-regression tests.
#
# Layout fingerprints are deterministic and platform-independent (pure
# geometry). Slow platform/font-dependent work (PNG rendering, md5 baselines,
# plotly conversion, and heavy device rendering) is opt-in and only runs when
# GGCHORD_RUN_SLOW_TESTS=1.

skip_unless_slow_tests <- function() {
  skip_if(Sys.getenv("GGCHORD_RUN_SLOW_TESTS") != "1",
          "set GGCHORD_RUN_SLOW_TESTS=1 to run slow rendering/integration tests")
}

#' Render a ggchord plot to a PNG file (fixed device), return the path
render_png <- function(p, file = tempfile(fileext = ".png"),
                       width = 8, height = 8) {
  grDevices::png(file, width = width, height = height,
                 units = "in", res = 150)
  on.exit(grDevices::dev.off())
  print(p)
  file
}

#' Deterministic fingerprint of a built ggchord layout (geometry + colors)
ggchord_layout_fingerprint <- function(p) {
  ggplot_build(p)
  l <- get_chord_layout()
  parts <- character(0)
  add_part <- function(name, value) {
    parts <<- c(parts, paste0(name, "=", value))
  }
  if (length(l$seq_arcs) > 0) {
    arcs <- do.call(rbind, l$seq_arcs)
    add_part("seq", paste(round(arcs$x, 6), round(arcs$y, 6),
                          arcs$seq_id, collapse = "|"))
  }
  if (!is.null(l$ribbon_polys)) {
    rp <- l$ribbon_polys
    add_part("ribbon", paste(round(rp$x, 6), round(rp$y, 6),
                             rp$group, rp$pident %||% "",
                             rp$fill %||% "", collapse = "|"))
  }
  if (nrow(l$gene_polys) > 0) {
    gp <- l$gene_polys
    add_part("gene", paste(round(gp$x, 6), round(gp$y, 6),
                           gp$group, collapse = "|"))
  }
  if (nrow(l$gene_labels) > 0) {
    gl <- l$gene_labels
    add_part("gene_label", paste(round(gl$text_x, 6), round(gl$text_y, 6),
                                 gl$text, collapse = "|"))
  }
  if (nrow(l$seq_labels_df) > 0) {
    sl <- l$seq_labels_df
    add_part("seq_label", paste(round(sl$text_x, 6), round(sl$text_y, 6),
                                sl$label, collapse = "|"))
  }
  if (nrow(l$axis_ticks) > 0) {
    at <- l$axis_ticks
    add_part("axis", paste(round(at$x0, 6), round(at$y0, 6),
                           round(at$x1, 6), round(at$y1, 6),
                           at$label, collapse = "|"))
  }
  add_part("seq_colors", paste(l$seq_colors, collapse = ","))
  add_part("ribbon_colors", paste(l$ribbon_colors, collapse = ","))
  add_part("gene_pal", paste(l$gene_pal, collapse = ","))
  paste(parts, collapse = ";;;")
}

#' The standard example plot used across visual-regression tests
example_ggchord_plot <- function() {
  data(seq_data_example, envir = environment())
  data(ribbon_data_example, envir = environment())
  data(gene_data_example, envir = environment())
  ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() + geom_ribbon() + geom_gene() + geom_axis()
}

#' Render the canonical set of regression plots to PNGs in a directory.
#' Returns a named list of plot name -> png path.
render_regression_set <- function(dir = tempfile()) {
  dir.create(dir, showWarnings = FALSE)
  data(seq_data_example, envir = environment())
  data(ribbon_data_example, envir = environment())
  data(gene_data_example, envir = environment())
  plots <- list(
    default = example_ggchord_plot(),
    seq_only = ggchord(seq_data_example) + geom_seq(),
    ribbon_pident = ggchord(seq_data_example, ribbon_data_example) +
      geom_seq() + geom_ribbon(),
    ribbon_query = ggchord(seq_data_example, ribbon_data_example) +
      geom_seq() + geom_ribbon(ribbon_color_scheme = "query"),
    ribbon_subject = ggchord(seq_data_example, ribbon_data_example) +
      geom_seq() + geom_ribbon(ribbon_color_scheme = "subject"),
    ribbon_single = ggchord(seq_data_example, ribbon_data_example) +
      geom_seq() +
      geom_ribbon(ribbon_color_scheme = "single", ribbon_colors = "orange"),
    gene_labels = ggchord(seq_data_example, gene_data = gene_data_example) +
      geom_seq() + geom_gene() + geom_gene_label(),
    gene_repel = ggchord(seq_data_example, gene_data = gene_data_example) +
      geom_seq() + geom_gene() + geom_gene_label_repel(seed = 123),
    axis_seq_labels = ggchord(seq_data_example) +
      geom_seq() + geom_axis() + geom_seq_label(),
    two_seq = {
      seq2 <- seq_data_example[seq_data_example$seq_id %in%
                                 c("MT108731.1", "MT118296.1"), ]
      rib2 <- ribbon_data_example[
        ribbon_data_example$qaccver %in% seq2$seq_id &
          ribbon_data_example$saccver %in% seq2$seq_id, ]
      gen2 <- gene_data_example[gene_data_example$seq_id %in% seq2$seq_id, ]
      ggchord(seq2, rib2, gen2) + geom_seq() + geom_ribbon() +
        geom_gene() + geom_axis()
    },
    three_seq = {
      seq3 <- seq_data_example[seq_data_example$seq_id %in%
                                 c("MT108731.1", "MT118296.1", "OQ646790.1"), ]
      rib3 <- ribbon_data_example[
        ribbon_data_example$qaccver %in% seq3$seq_id &
          ribbon_data_example$saccver %in% seq3$seq_id, ]
      gen3 <- gene_data_example[gene_data_example$seq_id %in% seq3$seq_id, ]
      ggchord(seq3, rib3, gen3) + geom_seq() + geom_ribbon() +
        geom_gene() + geom_axis()
    },
    theme_scale = ggchord(seq_data_example, ribbon_data_example,
                          gene_data_example) +
      geom_seq() + geom_ribbon() + geom_gene() + geom_axis() +
      theme(panel.background = element_rect(fill = "grey95")) +
      scale_color_manual(values = c("MT108731.1" = "#E41A1C",
                                    "MT118296.1" = "#377EB8",
                                    "OQ646790.1" = "#4DAF4A",
                                    "OR222515.1" = "#984EA3"))
  )
  out <- list()
  for (nm in names(plots)) {
    set.seed(123)
    out[[nm]] <- render_png(plots[[nm]], file = file.path(dir, paste0(nm, ".png")))
  }
  out
}
