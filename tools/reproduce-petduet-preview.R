# Reproduce the reported pETDuet-1 map and a SnapGene-reference configuration.
# Run from the repository root:
#   Rscript tools/reproduce-petduet-preview.R [output-directory]
# Uses the default view_ggchord() size; this is a focused visual acceptance
# script, not a package release check.

devtools::load_all(quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[1L] else
  file.path(tempdir(), "ggchord-petduet-preview")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path("data", "plasmid_example_pETDuet_1.rda"))
sequence <- plasmid_example_pETDuet_1
features <- find_common_features(sequence)
primers <- find_primer_bindings(sequence, set = "reference")
all_sites <- find_restriction_sites(sequence)
tracks <- position_feature_stack(spacing = .085,
  base_position = position_plasmid())

make_plot <- function(sites, boxed_primers) {
  primer_layers <- if (boxed_primers) {
    primer_labels <- transform(primers, label = name,
      feature_label_mode = "external", feature_label_fill = "#ECECEC",
      feature_label_colour = "#444444")
    list(geom_primer(data = primers),
      geom_feature_label_repel(data = primer_labels,
        label_layout = "callout", label_fit = "none",
        label_max_lines = 1L, max_overlaps = 0,
        position = position_plasmid(0)))
  } else list(geom_primer(data = primers),
    geom_primer_label_repel(data = primers, max_overlaps = 0))
  ggchord(sequence, validate = "none") +
    geom_seq(seq_style = "double") +
    geom_feature_plasmid(data = features, position = tracks) +
    geom_feature_label_repel(data = features, position = tracks,
      external = TRUE, max_overlaps = 0) +
    geom_restriction_site(data = sites) + primer_layers +
    geom_seq_center_label() + coord_circular(rotation = 90) +
    theme_ggchord_plasmid()
}

audit_plot <- function(plot) {
  # view_ggchord() has already frozen a device-sized layout on this plot.
  # Rebuilding via export_ggchord_layout() on R's null device would audit a
  # different physical label solution from the image just saved.
  layout <- plot$ggchord$ref$layout
  component <- function(name) {
    pieces <- lapply(names(layout$layer_geometry), function(id) {
      geometry <- layout$layer_geometry[[id]][[name]]
      if (!is.data.frame(geometry) || !nrow(geometry)) return(NULL)
      geometry$layer_id <- id
      geometry
    })
    pieces <- Filter(Negate(is.null), pieces)
    if (!length(pieces)) return(data.frame())
    ggchord:::ggchord_rbind_fill(pieces)
  }
  exported <- list(
    annotation_registry = layout$circular_annotation_registry,
    restriction = component("restriction_site"),
    labels = component("gene_label_repel")
  )
  all_boxes <- exported$annotation_registry
  all_boxes <- all_boxes[is.finite(all_boxes$bbox_xmin) &
    is.finite(all_boxes$bbox_xmax) &
    is.finite(all_boxes$bbox_ymin) &
    is.finite(all_boxes$bbox_ymax), , drop = FALSE]
  boxes <- all_boxes[all_boxes$side %in% "outer", , drop = FALSE]
  text_boxes <- data.frame(
    registry_id = all_boxes$registry_id, owner = all_boxes$owner,
    source_row = all_boxes$source_row,
    cx = (all_boxes$bbox_xmin + all_boxes$bbox_xmax) / 2,
    cy = (all_boxes$bbox_ymin + all_boxes$bbox_ymax) / 2,
    w = all_boxes$bbox_xmax - all_boxes$bbox_xmin,
    h = all_boxes$bbox_ymax - all_boxes$bbox_ymin,
    angle = 0
  )
  labels <- exported$labels
  if (nrow(labels) && all(c(".component", "text", "text_x", "text_y",
      "text_angle", "annotation_side") %in% names(labels))) {
    glyphs <- labels[labels$.component %in% "arc_text" &
      labels$annotation_side %in% "inner" &
      !is.na(labels$text) & nzchar(trimws(labels$text)) &
      is.finite(labels$text_x) & is.finite(labels$text_y) &
      is.finite(labels$text_angle), , drop = FALSE]
    if (nrow(glyphs)) {
      # The registry keeps one axis-aligned envelope for each whole label.
      # Around curved text its empty corners are not visible obstacles:
      # audit the final rendered glyph positions and rotations instead.
      measured <- ggchord:::ggchord_text_boxes(glyphs,
        units_per_inch = layout$text_units_per_inch, box_padding = 0)
      parent_key <- paste(all_boxes$owner, all_boxes$source_row, sep = "\r")
      glyph_key <- paste(glyphs$layer_id, glyphs$source_row, sep = "\r")
      parent <- match(glyph_key, parent_key)
      glyph_boxes <- data.frame(
        registry_id = ifelse(is.na(parent),
          paste0("arc-label-", glyphs$layer_id, "-", glyphs$source_row),
          all_boxes$registry_id[parent]),
        owner = glyphs$layer_id, source_row = glyphs$source_row,
        cx = measured$cx, cy = measured$cy, w = measured$w,
        h = measured$h, angle = measured$angle
      )
      replace <- all_boxes$side %in% "inner" &
        all_boxes$kind %in% "label" & parent_key %in% glyph_key
      text_boxes <- rbind(text_boxes[!replace, , drop = FALSE], glyph_boxes)
    }
  }
  overlap <- 0L
  if (nrow(boxes) > 1L) {
    pairs <- utils::combn(seq_len(nrow(boxes)), 2L)
    overlap <- sum(vapply(seq_len(ncol(pairs)), function(k) {
      i <- pairs[1L, k]; j <- pairs[2L, k]
      min(boxes$bbox_xmax[c(i, j)]) >
        max(boxes$bbox_xmin[c(i, j)]) + 1e-6 &&
        min(boxes$bbox_ymax[c(i, j)]) >
        max(boxes$bbox_ymin[c(i, j)]) + 1e-6
    }, logical(1L)))
  }
  site_labels <- exported$restriction
  site_labels <- site_labels[site_labels$.component %in% "label" &
    !is.na(site_labels$label), , drop = FALSE]
  hits_box <- function(segment, box) {
    # Intersect in the measured box's own coordinate system so a diagonal
    # glyph is not replaced by the larger axis-aligned rectangle around it.
    cos_a <- cos(box$angle)
    sin_a <- sin(box$angle)
    rel_x <- segment$x1 - box$cx
    rel_y <- segment$y1 - box$cy
    dx <- segment$x2 - segment$x1
    dy <- segment$y2 - segment$y1
    origin <- c(rel_x * cos_a + rel_y * sin_a,
      -rel_x * sin_a + rel_y * cos_a)
    delta <- c(dx * cos_a + dy * sin_a,
      -dx * sin_a + dy * cos_a)
    half <- c(box$w, box$h) / 2
    lower <- 1e-5
    upper <- 1 - 1e-5
    for (axis in seq_len(2L)) {
      if (abs(delta[axis]) <= 1e-9) {
        if (origin[axis] <= -half[axis] ||
            origin[axis] >= half[axis]) return(FALSE)
      } else {
        interval <- sort(c((-half[axis] - origin[axis]) / delta[axis],
          (half[axis] - origin[axis]) / delta[axis]))
        lower <- max(lower, interval[1L])
        upper <- min(upper, interval[2L])
      }
    }
    lower + 1e-7 < upper
  }
  leaders <- list()
  restriction <- exported$restriction
  restriction <- restriction[restriction$restriction_component %in% "leader" &
    is.finite(restriction$source_row), , drop = FALSE]
  if (nrow(restriction)) {
    key <- paste(restriction$layer_id, restriction$source_row)
    for (name in unique(key)) {
      path <- restriction[key == name, , drop = FALSE]
      if (nrow(path) < 2L) next
      for (k in seq_len(nrow(path) - 1L)) {
        leaders[[length(leaders) + 1L]] <- data.frame(
          owner = path$layer_id[k], source_row = path$source_row[k],
          x1 = path$x[k], y1 = path$y[k],
          x2 = path$x[k + 1L], y2 = path$y[k + 1L]
        )
      }
    }
  }
  if (nrow(labels) && all(c(".component", "xend", "yend") %in%
      names(labels))) {
    segments <- labels[labels$.component %in% "segment" &
      is.finite(labels$source_row), , drop = FALSE]
    for (k in seq_len(nrow(segments))) {
      leaders[[length(leaders) + 1L]] <- data.frame(
        owner = segments$layer_id[k],
        source_row = segments$source_row[k],
        x1 = segments$x[k], y1 = segments$y[k],
        x2 = segments$xend[k], y2 = segments$yend[k]
      )
    }
  }
  leaders <- Filter(function(segment) {
    coordinates <- unlist(segment[c("x1", "y1", "x2", "y2")],
      use.names = FALSE)
    all(is.finite(coordinates)) &&
      (segment$x2 - segment$x1)^2 +
        (segment$y2 - segment$y1)^2 > 1e-16
  }, leaders)
  source_key <- function(segment) {
    paste(segment$owner, segment$source_row, sep = "\r")
  }
  crosses <- function(first, second) {
    rx <- first$x2 - first$x1
    ry <- first$y2 - first$y1
    sx <- second$x2 - second$x1
    sy <- second$y2 - second$y1
    denominator <- rx * sy - ry * sx
    if (abs(denominator) <= 1e-12) return(FALSE)
    dx <- second$x1 - first$x1
    dy <- second$y1 - first$y1
    t <- (dx * sy - dy * sx) / denominator
    u <- (dx * ry - dy * rx) / denominator
    # Shared endpoints and a route's own bend are not line crossings.
    t > 1e-9 && t < 1 - 1e-9 && u > 1e-9 && u < 1 - 1e-9
  }
  crossing_pairs <- list()
  if (length(leaders) > 1L) {
    for (i in seq_len(length(leaders) - 1L)) {
      first <- leaders[[i]]
      for (j in seq.int(i + 1L, length(leaders))) {
        second <- leaders[[j]]
        if (identical(source_key(first), source_key(second)) ||
            !crosses(first, second)) next
        sources <- sort(c(source_key(first), source_key(second)))
        key <- paste(sources, collapse = "\n")
        crossing_pairs[[key]] <- paste(gsub("\r", " ", sources,
          fixed = TRUE), collapse = " <-> ")
      }
    }
  }
  leader_text_hits <- character()
  for (leader in leaders) {
    for (k in seq_len(nrow(text_boxes))) {
      box <- text_boxes[k, , drop = FALSE]
      if (identical(as.character(leader$owner),
          as.character(box$owner)) &&
          identical(as.integer(leader$source_row),
            as.integer(box$source_row))) next
      if (hits_box(leader, box)) {
        leader_text_hits <- c(leader_text_hits,
          paste(leader$owner, leader$source_row, box$registry_id))
      }
    }
  }
  if (length(crossing_pairs)) print(unname(unlist(crossing_pairs)))
  if (length(leader_text_hits)) {
    print(unique(leader_text_hits))
  }
  data.frame(
    sites = nrow(site_labels), outer_labels = nrow(boxes),
    bbox_overlaps = overlap,
    leader_crossings = length(crossing_pairs),
    leader_text_hits = length(unique(leader_text_hits))
  )
}

configurations <- list(
  reported = list(
    sites = filter_restriction_sites(all_sites, set = "unique_6plus",
      parent_set = "commercial_nonredundant"),
    boxed_primers = TRUE, expected_sites = 54L),
  reference_enzymes = list(
    sites = filter_restriction_sites(all_sites, set = "unique_6plus",
      parent_set = "commercial_all"),
    boxed_primers = FALSE, expected_sites = 58L)
)
results <- list()
for (name in names(configurations)) {
  config <- configurations[[name]]
  plot <- make_plot(config$sites, config$boxed_primers)
  image_path <- view_ggchord(plot, viewer = "none")
  destination <- file.path(output_dir, paste0("pETDuet-1-", name,
    ".png"))
  stopifnot(file.copy(image_path, destination, overwrite = TRUE))
  result <- audit_plot(plot)
  result$configuration <- name
  result$image <- destination
  print(result, row.names = FALSE)
  stopifnot(result$sites == config$expected_sites,
    result$bbox_overlaps == 0L, result$leader_crossings == 0L,
    result$leader_text_hits == 0L)
  results[[name]] <- result
}
print(do.call(rbind, results), row.names = FALSE)
