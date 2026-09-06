# Reproducible release acceptance outside the fast testthat suite.
# Rscript tools/validate-release.R geometry|links|labels|output|benchmark|examples [directory]
# Outputs are temporary by default; no README/site figures are overwritten.
suppressPackageStartupMessages(library(ggplot2))
pkgload::load_all(quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "geometry"
out <- if (length(args) > 1) args[2] else tempfile("ggchord-acceptance-")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
message("Output: ", normalizePath(out))
writeLines(capture.output(sessionInfo()), file.path(out, "session.txt"))

with_device <- function(code, width = 11, height = 7) {
  grDevices::pdf(NULL, width = width, height = height)
  on.exit(grDevices::dev.off())
  force(code)
}
finite_geometry <- function(x) {
  if (is.data.frame(x)) {
    cols <- intersect(c("x", "y", "x0", "y0", "x1", "y1", "text_x", "text_y"), names(x))
    stopifnot(all(vapply(x[cols], function(z) all(is.finite(z)), logical(1))))
  } else if (is.list(x)) invisible(lapply(x, finite_geometry))
}
label_metrics <- function(layout, expected) {
  labels <- layout$gene_labels
  stopifnot(sum(!is.na(labels$text) & nzchar(labels$text)) == expected)
  stopifnot(!ggchord:::ggchord_label_box_conflicts(labels,
    units_per_inch = layout$text_units_per_inch, box_padding = 0))
  paths <- layout$gene_label_segments
  if (nrow(paths) > 1L) for (i in seq_len(nrow(paths) - 1L)) {
    j <- seq.int(i + 1L, nrow(paths))
    j <- j[paths$group[j] != paths$group[i]]
    stopifnot(!any(ggchord:::ggchord_segments_cross(paths$x0[i], paths$y0[i],
      paths$x1[i], paths$y1[i], paths$x0[j], paths$y0[j], paths$x1[j], paths$y1[j])))
  }
  if (identical(layout$gene_label_layout, "auto")) {
    frame <- ggchord:::ggchord_label_curve_frame(labels, layout$seq_arcs)
    columns <- abs(frame$outward_x) >= abs(frame$outward_y)
    side <- ifelse(frame$outward_x < 0, "left", "right")
    stopifnot(any(columns), any(!columns))
    for (direction in unique(side[columns])) {
      stopifnot(length(unique(round(labels$text_x[columns & side == direction], 8))) == 1)
    }
    stopifnot(all(is.na(labels$.radial_parameter[columns])),
      all(is.finite(labels$.radial_parameter[!columns])))
  }
  data.frame(labels = expected, conflicts = 0, crossings = 0,
    max_track = max(labels$label_track, na.rm = TRUE))
}
base_plot <- function() ggchord(seq_data_example, ribbon_data_example,
  gene_data_example, validate = "none") + geom_seq() + geom_link_ribbon() + geom_gene()

if (mode == "geometry") {
  sequences <- data.frame(accver = c("A", "B"), length = 1000)
  features <- data.frame(accver = rep(c("A", "B"), each = 4),
    start = rep(c(100, 300, 500, 700), 2), end = rep(c(160, 360, 560, 760), 2),
    strand = rep(c("+", "-"), 4), type = rep(c("arrow", "block", "chevron", "lollipop"), 2))
  genes <- transform(features, anno = type)
  ribbons <- data.frame(qaccver = "A", saccver = "B", qstart = 100, qend = 160,
    sstart = 500, send = 560, length = 61, pident = 95)
  results <- list()
  for (curvature in c(0.5, 0.8, 1, 1.2)) {
    p <- ggchord(sequences, ribbons, genes, validate = "none") +
      geom_seq(seq_radius = c(2.5, 1.8), seq_curvature = c(curvature, 1),
        seq_orientation = c(1, -1), seq_gap = c(0.04, 0.08)) +
      geom_gene(gene_offset = 0.15) +
      geom_feature(aes(feature_shape = type), data = features, feature_offset = -0.15) +
      scale_feature_shape_manual(values = setNames(unique(features$type), unique(features$type))) +
      geom_seq_region(data = features[1, ], region_side = "outside") +
      geom_link_ribbon() +
      geom_gene_label(size = 2, gene_label_overlap = "allow") +
      coord_chord(rotation = 35)
    layout <- with_device(get_chord_layout(p))
    finite_geometry(c(layout$seq_arcs, list(layout$gene_polys, layout$ribbon_polys,
      layout$region_polys, layout$axis_lines, layout$axis_ticks)))
    exported <- with_device(export_ggchord_layout(p))
    stopifnot(nrow(exported$feature) > 0, nrow(exported$gene) > 0,
      nrow(exported$axis) > 0, nrow(exported$ribbon) > 0)
    stopifnot(length(unique(exported$feature$source_row[exported$feature$component == "gene_poly"])) == 8)
    finite_geometry(exported$feature)
    ggplot2::ggsave(file.path(out, paste0("curvature-", curvature, ".png")), p,
      width = 11, height = 7, dpi = 100)
    results[[length(results) + 1]] <- data.frame(curvature, finite = TRUE, features = 8)
  }
  write.csv(do.call(rbind, results), file.path(out, "geometry.csv"), row.names = FALSE)
} else if (mode == "links") {
  sequences <- data.frame(accver = c("A", "B", "C"), length = 1000)
  links <- data.frame(qaccver = "A", saccver = c("B", "C"),
    qpos = 200, spos = 500)
  ribbons <- transform(links, qstart = 100, qend = 400,
    sstart = 600, send = 400)
  features <- data.frame(accver = "A", start = c(100, 300, 500, 700),
    end = c(200, 400, 600, 800), strand = c("-", "-", "+", "+"),
    type = c("CDS", "tRNA", "repeat", "promoter"))
  shapes <- c(CDS = "arrow", tRNA = "block", "repeat" = "chevron",
    promoter = "lollipop")
  base <- ggchord(sequences) + geom_seq()
  plots <- list(
    line_query = base + geom_link_line(data = links, link_branch = "query",
      arrow = arrow(length = unit(2, "mm"))),
    line_subject = base + geom_link_line(
      data = transform(links, qaccver = saccver, saccver = qaccver,
        qpos = spos, spos = qpos),
      link_branch = "subject", arrow = arrow(ends = "both", length = unit(2, "mm"))),
    ribbon_query = ggchord(sequences, ribbons, validate = "none") + geom_seq() +
      geom_link_ribbon(fill = "#3B90B6", link_branch = "query"),
    feature_keys = base +
      geom_feature(aes(feature_shape = type), data = features) +
      scale_feature_shape_manual(values = shapes),
    smooth = ggchord(sequences, ribbons, validate = "none") + geom_seq() +
      geom_feature(aes(feature_shape = type), data = features) +
      scale_feature_shape_manual(values = shapes) +
      geom_link_ribbon(fill = "#3B90B6", link_avoid = "smooth"),
    uniform = ggchord(sequences, ribbons, validate = "none") + geom_seq() +
      geom_feature(aes(feature_shape = type), data = features) +
      scale_feature_shape_manual(values = shapes) +
      geom_link_ribbon(fill = "#3B90B6", link_avoid = "uniform")
  )
  for (name in names(plots)) {
    layout <- with_device(get_chord_layout(plots[[name]]), 8, 6)
    finite_geometry(layout)
    for (extension in c("png", "pdf", "svg")) {
      ggplot2::ggsave(file.path(out, paste0(name, ".", extension)),
        plots[[name]], width = 8, height = 6, dpi = 120)
    }
  }
} else if (mode == "labels") {
  cases <- list(
    default = list(radius = rep(2.5, 4), curvature = rep(1, 4), orientation = rep(1, 4), rotation = 45),
    unequal = list(radius = c(3.3, 2.5, 1.8, 1.25), curvature = c(.8, 1.2, .7, 1.1), orientation = c(1,-1,1,-1), rotation = 35))
  results <- list()
  for (layout_mode in c("radial", "auto")) for (name in names(cases)) {
    x <- cases[[name]]
    p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example, validate = "none") +
      geom_seq(seq_radius = x$radius, seq_curvature = x$curvature, seq_orientation = x$orientation) +
      geom_link_ribbon() + geom_gene() + geom_gene_label_repel(gene_label_layout = layout_mode) +
      coord_chord(rotation = x$rotation)
    elapsed <- system.time(layout <- with_device(get_chord_layout(p)))[["elapsed"]]
    metrics <- label_metrics(layout, nrow(gene_data_example))
    results[[length(results) + 1L]] <- cbind(mode = layout_mode, case = name, elapsed, metrics)
    saveRDS(layout, file.path(out, paste0(layout_mode, "-", name, ".rds")))
    ggplot2::ggsave(file.path(out, paste0(layout_mode, "-", name, ".png")), p,
      width = 11, height = 7, dpi = 100)
    message(layout_mode, " / ", name, ": ", elapsed, "s; no conflicts or crossings")
  }
  write.csv(do.call(rbind, results), file.path(out, "labels.csv"), row.names = FALSE)
} else if (mode == "output") {
  stopifnot(requireNamespace("png", quietly = TRUE), requireNamespace("svglite", quietly = TRUE))
  p <- base_plot() + geom_gene_label_repel()
  original_device <- grDevices::dev.cur()
  for (device in c("png", "svg")) {
    file <- view_ggchord(p, device = device, viewer = "none")
    file.copy(file, file.path(out, paste0("default.", device)))
    stopifnot(grDevices::dev.cur() == original_device)
  }
  for (units in c("in", "cm", "mm", "px")) {
    multiplier <- c("in" = 1, cm = 2.54, mm = 25.4, px = 100)[units]
    file <- view_ggchord(p, width = 11 * multiplier, height = 7 * multiplier,
      units = units, dpi = 100, viewer = "none")
    stopifnot(identical(dim(png::readPNG(file))[1:2], c(700L, 1100L)),
      grDevices::dev.cur() == original_device)
    file.copy(file, file.path(out, paste0("explicit-", units, ".png")))
  }
  for (position in c("none", "top", "bottom", "left", "right")) {
    sample <- ggchord(seq_data_example, ribbon_data_example, validate = "none") +
      geom_seq() + geom_link_ribbon() +
      scale_seq_colour_manual(values = setNames(c("red", "blue", "green", "orange"), seq_data_example$accver),
        labels = paste("A longer sequence description", seq_len(4))) +
      theme(legend.position = position)
    if (position %in% c("top", "bottom")) {
      sample <- sample + guides(seq_colour = guide_ggchord_legend(nrow = 2))
    }
    file <- view_ggchord(sample, width = 6, viewer = "none")
    file.copy(file, file.path(out, paste0("legend-", position, ".png")))
  }
  ggplot2::ggsave(file.path(out, "explicit.pdf"), p, width = 11, height = 7)
  stopifnot(grDevices::dev.cur() == original_device)
} else if (mode == "benchmark") {
  results <- list()
  for (n in c(32, 80)) {
    sequences <- data.frame(accver = "A", length = 100000)
    starts <- seq(100, 99000, length.out = n)
    genes <- data.frame(accver = "A", start = starts, end = starts + 50,
      strand = "+", anno = paste0("gene", seq_len(n)))
    p <- ggchord(sequences, gene_data = genes, validate = "none") + geom_seq() +
      geom_gene() + geom_gene_label_repel(size = 2)
    elapsed <- system.time(layout <- with_device(get_chord_layout(p), 12, 12))[["elapsed"]]
    results[[length(results) + 1]] <- cbind(n, elapsed, label_metrics(layout, n))
    message(n, " labels: ", elapsed, "s")
  }
  write.csv(do.call(rbind, results), file.path(out, "benchmark.csv"), row.names = FALSE)
} else if (mode == "examples") {
  for (file in list.files("vignettes", "\\.Rmd$", full.names = TRUE)) {
    lines <- readLines(file)
    starts <- grep("^```\\{r", lines)
    env <- new.env(parent = globalenv())
    for (start in starts) {
      header <- lines[start]
      # Image inclusion is checked by R CMD build; import paths are repository-only.
      if (grepl("-img|setup|import-helpers", header)) next
      end <- start + which(lines[-seq_len(start)] == "```")[1]
      code <- parse(text = lines[seq.int(start + 1, end - 1)])
      message(basename(file), ": ", header)
      with_device(for (expression in code) {
        value <- eval(expression, env)
        if (inherits(value, "ggplot")) invisible(ggplot2::ggplotGrob(value))
      })
    }
  }
} else stop("Unknown validation mode: ", mode)
message("PASS: ", mode)
