# Explicit local fitting checks; no testthat or package release checks.
pkgload::load_all(quiet = TRUE)
local({
  previous <- grDevices::dev.cur()
  grDevices::pdf(NULL, width = 7, height = 5)
  owned <- grDevices::dev.cur()
  on.exit({
    if (owned %in% grDevices::dev.list()) grDevices::dev.off(owned)
    if (previous %in% grDevices::dev.list()) grDevices::dev.set(previous)
  }, add = TRUE)
  sequences <- data.frame(accver = "A", length = 1000)
  features <- data.frame(accver = "A", start = 100, end = 350,
    strand = "+", type = "CDS", label = "LONG GENE LABEL")
  make_plot <- function(coordinate) ggchord(sequences, validate = "none") +
    geom_seq() + geom_feature(data = features, width = .12) +
    geom_feature_label_repel(data = features, label_fit = "shrink",
      label_min_size = 1.5, size = 5) + coordinate
  for (coordinate in list(coord_collinear(), coord_circular(), coord_chord())) {
    exported <- export_ggchord_layout(make_plot(coordinate), include = "labels")
    geometry <- exported$labels
    text <- geometry[geometry$.component %in% "text", , drop = FALSE]
    stopifnot(nrow(text) > 0L, all(text$label_text_original == features$label),
      all(text$size >= 1.5 & text$size <= 5), any(text$size < 5),
      all(text$label_size_final == text$size),
      nrow(exported$label_fit_diagnostics) == nrow(features),
      setequal(exported$label_fit_diagnostics$layer_id, unique(text$layer_id)))
  }
  # A rendered text grob must retain its resolved size rather than have the
  # original fixed size reapplied by the composite text-style wrapper.
  collect <- function(grob) {
    if (inherits(grob, "text") && any(as.character(grob$label) == features$label))
      return(list(grob))
    children <- if (inherits(grob, "gtable")) grob$grobs else grob$children
    unlist(lapply(children, collect), recursive = FALSE)
  }
  text <- collect(ggplot2::ggplotGrob(make_plot(coord_collinear())))
  stopifnot(length(text) > 0L, all(vapply(text, function(grob)
    all(grob$gp$fontsize >= 1.5 * 72.27 / 25.4 & grob$gp$fontsize < 5 * 72.27 / 25.4), logical(1))))

  short <- transform(features, end = 102,
    label = "A long annotation that cannot fit inside a two-base feature")
  fallback <- ggchord(sequences, validate = "none") + geom_seq() +
    geom_feature(data = short, width = .08) +
    geom_feature_label_repel(data = short, label_fit = "shrink", size = 5,
      external = FALSE) + coord_circular()
  layout <- get_chord_layout(fallback)
  fit <- ggchord_rbind_fill(lapply(layout$layer_layouts, `[[`, "label_fit_diagnostics"))
  fit <- fit[fit$label_text_original == short$label, , drop = FALSE]
  stopifnot(nrow(fit) > 0L, all(fit$label_size_final >= 1.5),
    any(grepl("minimum_size_exceeded", fit$label_fit_reason)))

  genes <- transform(features, anno = label)
  gene_plot <- ggchord(sequences, gene_data = genes, validate = "none") +
    geom_seq() + geom_gene(width = .12) +
    geom_gene_label_repel(label_fit = "shrink", size = 5) + coord_collinear()
  text <- export_ggchord_layout(gene_plot, include = "labels")$labels
  text <- text[text$.component %in% "text", , drop = FALSE]
  stopifnot(any(text$size < 5), all(text$label_text_original == genes$anno))
})
cat("Explicit bounded label shrinking, rendered fonts and source text passed.\n")
