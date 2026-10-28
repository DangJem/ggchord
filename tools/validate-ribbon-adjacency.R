# Direct adjacency-style checks and optional PNG/PDF/SVG acceptance.
# Rscript tools/validate-ribbon-adjacency.R [output-directory]
# Four sequences distinguish a circular wrap neighbour from a skipped track.
pkgload::load_all(quiet = TRUE)
seqs <- data.frame(accver = c("A", "B", "C", "D"), length = 1000)
links <- data.frame(qaccver = c("A", "B", "C", "A", "A", "B"),
  saccver = c("B", "C", "D", "C", "D", "D"),
  qstart = c(100, 350, 600, 350, 650, 650),
  qend = c(250, 500, 750, 500, 800, 800),
  sstart = c(100, 350, 600, 100, 100, 350),
  send = c(250, 500, 750, 250, 250, 500))
base <- ggchord(seqs, validate = "none")
make_plot <- function(coordinate, order = seqs$accver, adjacent_only = NULL) {
  base + geom_link_ribbon(data = links, adjacent_only = adjacent_only,
    fill = "#91AFC0", colour = "#3E4550", linewidth = .7, alpha = .4) +
    geom_seq(seq_order = order) + geom_seq_label() +
    scale_seq_position_continuous(breaks = NULL, minor_breaks = NULL) +
    coordinate
}
plots <- list(
  linear_default = make_plot(coord_collinear()),
  linear_all = make_plot(coord_collinear(), adjacent_only = FALSE),
  linear_reordered = make_plot(coord_collinear(), c("A", "C", "B", "D"), FALSE),
  chord_wrap = make_plot(coord_chord())
)
# Every branch exists before any build. Direct export must still resolve its
# own layout rather than reuse a sibling branch's cached source-layer rows.
for (name in c(names(plots), rev(names(plots)))) {
  exported <- export_ggchord_layout(plots[[name]], include = "ribbon")$ribbon
  expected <- if (name == "linear_default") 1:3 else 1:6
  stopifnot(setequal(unique(exported$source_row), expected))
}
styles <- function(plot) {
  # Source identities belong to the layout export; ggplot2 discards those
  # non-aesthetic columns. Build and associate the rendered groups with the
  # source identities from that same layout.
  drawn <- ggplot2::ggplot_build(plot)$data[[1L]]
  exported <- export_ggchord_layout(plot, include = "ribbon")$ribbon
  sources <- unique(exported[c("group", "source_row")])
  rendered <- unique(drawn[c("group", "linetype")])
  result <- merge(sources, rendered, by = "group", all.x = TRUE)
  stopifnot(!anyNA(result$linetype))
  unique(result[c("source_row", "linetype")])
}
assert_styles <- function(plot, rows, expected) {
  actual <- styles(plot)
  stopifnot(setequal(actual$source_row, rows),
    identical(as.character(actual$linetype[match(rows, actual$source_row)]), expected))
}
assert_styles(plots$linear_default, 1:3, rep("solid", 3))
assert_styles(plots$linear_all, 1:6,
  c("solid", "solid", "solid", "dashed", "dashed", "dashed"))
assert_styles(plots$linear_reordered, 1:6,
  c("dashed", "solid", "dashed", "solid", "dashed", "solid"))
assert_styles(plots$chord_wrap, 1:6,
  c("solid", "solid", "solid", "dashed", "solid", "dashed"))

# Explicit layer styles and mapped values with user scales own their line type.
fixed <- base + geom_link_ribbon(data = links, adjacent_only = FALSE,
  linetype = "dotted") + geom_seq() + coord_collinear()
assert_styles(fixed, 1:6, rep("dotted", 6))
mapped <- base + geom_link_ribbon(ggplot2::aes(linetype = "chosen"),
  data = links, adjacent_only = FALSE) + geom_seq() + coord_collinear() +
  ggplot2::scale_linetype_manual(values = c(chosen = "dotdash"))
assert_styles(mapped, 1:6, rep("dotdash", 6))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) && nzchar(args[1L])) {
  out <- args[1L]
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  original <- grDevices::dev.cur()
  for (name in names(plots)) for (extension in c("png", "pdf", "svg")) {
    file <- file.path(out, paste0(name, ".", extension))
    ggplot2::ggsave(file, plots[[name]], width = 9, height = 7, dpi = 140,
      device = if (extension == "pdf") grDevices::pdf else extension)
    stopifnot(file.exists(file), file.info(file)$size > 0,
      identical(grDevices::dev.cur(), original))
  }
  writeLines(capture.output(sessionInfo()), file.path(out, "session.txt"))
}
cat("Ribbon adjacency, display order, circular wrap and user styles passed.\n")
