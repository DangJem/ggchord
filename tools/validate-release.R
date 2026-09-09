# Reproducible release acceptance outside the fast testthat suite.
# Rscript tools/validate-release.R genome|feature|geometry|links|labels|output|benchmark|examples [directory]
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
feature_polygons_overlap <- function(a, b) {
  if (max(a$x) < min(b$x) || max(b$x) < min(a$x) ||
      max(a$y) < min(b$y) || max(b$y) < min(a$y)) return(FALSE)
  inside <- any(vapply(seq_len(nrow(a)), function(i) {
    ggchord:::ggchord_point_in_polygon(a$x[i], a$y[i], b)
  }, logical(1))) || any(vapply(seq_len(nrow(b)), function(i) {
    ggchord:::ggchord_point_in_polygon(b$x[i], b$y[i], a)
  }, logical(1)))
  if (inside) return(TRUE)
  aa <- rbind(cbind(a$x, a$y), c(a$x[1], a$y[1]))
  bb <- rbind(cbind(b$x, b$y), c(b$x[1], b$y[1]))
  for (i in seq_len(nrow(aa) - 1L)) {
    for (j in seq_len(nrow(bb) - 1L)) {
      if (ggchord:::ggchord_segment_intersects(
          aa[i, ], aa[i + 1L, ], bb[j, ], bb[j + 1L, ])) return(TRUE)
    }
  }
  FALSE
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

if (mode == "genome") {
  p <- ggchord(single_genome_example, gene_data = single_gene_example,
               validate = "none") +
    geom_seq(seq_style = "double") + geom_gene(position = "plasmid") +
    geom_gene_label_repel(position = "plasmid") +
    geom_restriction_site(data = restriction_site_example) +
    geom_seq_center_label() + coord_circular(gap = 10)
  layout <- with_device(get_chord_layout(p), 8, 7)
  finite_geometry(list(layout$seq_arcs, layout$gene_polys,
                       layout$restriction_sites))
  stopifnot(nrow(layout$restriction_sites) > nrow(restriction_site_example))
  exported <- with_device(export_ggchord_layout(
    p, include = c("seq", "gene", "restriction"), original_data = TRUE
  ), 8, 7)
  stopifnot(exported$metadata$coordinate == "circular",
            exported$metadata$gap == 10,
            length(unique(stats::na.omit(exported$restriction$source_row))) ==
              nrow(restriction_site_example))
  ggplot2::ggsave(file.path(out, "single-genome.png"), p,
                  width = 8, height = 7, dpi = 120)
  ggplot2::ggsave(file.path(out, "single-genome.pdf"), p,
                  width = 8, height = 7)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(file.path(out, "single-genome.svg"), p,
                    width = 8, height = 7)
  }
} else if (mode == "feature") {
  objects <- c("plasmid_example_pUC19c",
    "plasmid_example_pBluescript_II_SK_plus")
  files <- c("pUC19c-feature", "pBluescript-feature")
  for (i in seq_along(objects)) {
    utils::data(list = objects[i])
    sequence <- get(objects[i])
    features <- find_common_features(sequence)
    tracks <- position_feature_stack(
      spacing = .10, base_position = position_plasmid()
    )
    p <- ggchord(sequence) + geom_seq() +
      geom_feature_plasmid(data = features, position = tracks) +
      geom_feature_label_repel(data = features, position = tracks) +
      geom_seq_center_label() + scale_seq_position_continuous() +
      coord_circular(rotation = 90) + theme_ggchord_plasmid()
    layout <- with_device(get_chord_layout(p), 7, 7)
    finite_geometry(list(layout$gene_polys, layout$gene_labels,
      layout$gene_label_segments))
    stopifnot(all(layout$gene_labels$feature_label_mode %in%
      c("inside", "adjacent", "callout")),
      any(layout$gene_labels$feature_label_mode == "inside"),
      any(layout$gene_labels$feature_label_mode != "inside"))
    stopifnot(all(ggchord:::ggchord_label_conflict_counts(
      layout$gene_labels, units_per_inch = layout$text_units_per_inch,
      box_padding = .005
    ) == 0L))
    boxes <- ggchord:::ggchord_text_boxes(
      layout$gene_labels, units_per_inch = layout$text_units_per_inch,
      box_padding = .04
    )
    polygon_rows <- layout$gene_polys$.component == "polygon"
    polygons <- split(layout$gene_polys[polygon_rows, , drop = FALSE],
      layout$gene_polys$group[polygon_rows])
    for (j in seq_len(nrow(layout$gene_labels))) {
      stopifnot(!ggchord:::ggchord_label_hits_features(
        boxes[j, , drop = FALSE], polygons,
        source_row = layout$gene_labels$source_row[j],
        allow_own = layout$gene_labels$feature_label_mode[j] == "inside"
      ))
    }
    if (length(polygons) > 1L) {
      for (j in seq_len(length(polygons) - 1L)) {
        for (k in seq.int(j + 1L, length(polygons))) {
          if (unique(polygons[[j]]$source_row)[1L] !=
              unique(polygons[[k]]$source_row)[1L]) {
            source_j <- unique(polygons[[j]]$source_row)[1L]
            source_k <- unique(polygons[[k]]$source_row)[1L]
            interval_overlap <- max(
              features$start[source_j], features$start[source_k]
            ) <= min(features$end[source_j], features$end[source_k])
            if (interval_overlap &&
                feature_polygons_overlap(polygons[[j]], polygons[[k]])) {
              stop("feature polygons overlap: ",
                source_j, " / ", source_k)
            }
          }
        }
      }
    }
    stopifnot(all(layout$gene_labels$text_angle <= 90 |
      layout$gene_labels$text_angle >= 270))
    ggplot2::ggsave(file.path(out, paste0(files[i], ".png")), p,
      width = 7, height = 7, dpi = 150)
    ggplot2::ggsave(file.path(out, paste0(files[i], ".pdf")), p,
      width = 7, height = 7)
    if (requireNamespace("svglite", quietly = TRUE)) {
      ggplot2::ggsave(file.path(out, paste0(files[i], ".svg")), p,
        width = 7, height = 7)
    }
  }
} else if (mode == "geometry") {
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
