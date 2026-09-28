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
  exported <- export_ggchord_layout(plot,
    include = c("feature", "labels", "restriction"))
  all_boxes <- exported$annotation_registry
  all_boxes <- all_boxes[is.finite(all_boxes$bbox_xmin) &
    is.finite(all_boxes$bbox_xmax) &
    is.finite(all_boxes$bbox_ymin) &
    is.finite(all_boxes$bbox_ymax), , drop = FALSE]
  boxes <- all_boxes[all_boxes$side == "outer", , drop = FALSE]
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
  site_labels <- site_labels[site_labels$.component == "label" &
    !is.na(site_labels$label), , drop = FALSE]
  hits_box <- function(segment, box) {
    lower <- 1e-5
    upper <- 1 - 1e-5
    for (axis in c("x", "y")) {
      start <- segment[[paste0(axis, "1")]]
      delta <- segment[[paste0(axis, "2")]] - start
      edge <- box[[paste0("bbox_", axis, "min")]]
      limit <- box[[paste0("bbox_", axis, "max")]]
      if (abs(delta) <= 1e-9) {
        if (start <= edge || start >= limit) return(FALSE)
      } else {
        interval <- sort(c((edge - start) / delta,
          (limit - start) / delta))
        lower <- max(lower, interval[1L])
        upper <- min(upper, interval[2L])
      }
    }
    lower + 1e-7 < upper
  }
  leaders <- list()
  restriction <- exported$restriction
  restriction <- restriction[restriction$restriction_component == "leader" &
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
  labels <- exported$labels
  if (nrow(labels) && all(c(".component", "xend", "yend") %in%
      names(labels))) {
    segments <- labels[labels$.component == "segment" &
      labels$annotation_side %in% "outer" &
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
  leader_text_hits <- character()
  for (leader in leaders) {
    for (k in seq_len(nrow(all_boxes))) {
      box <- all_boxes[k, , drop = FALSE]
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
  data.frame(
    sites = nrow(site_labels), outer_labels = nrow(boxes),
    bbox_overlaps = overlap,
    leader_crossings = sum(boxes$leader_crossing_count, na.rm = TRUE) / 2,
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
  stopifnot(result$sites == config$expected_sites,
    result$bbox_overlaps == 0L, result$leader_crossings == 0L,
    result$leader_text_hits == 0L)
  results[[name]] <- result
}
print(do.call(rbind, results), row.names = FALSE)
