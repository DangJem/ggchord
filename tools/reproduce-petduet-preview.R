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
  boxes <- exported$annotation_registry
  boxes <- boxes[boxes$side == "outer" &
    is.finite(boxes$bbox_xmin) & is.finite(boxes$bbox_xmax) &
    is.finite(boxes$bbox_ymin) & is.finite(boxes$bbox_ymax), ,
    drop = FALSE]
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
  data.frame(
    sites = nrow(site_labels), outer_labels = nrow(boxes),
    bbox_overlaps = overlap,
    leader_crossings = sum(boxes$leader_crossing_count, na.rm = TRUE) / 2
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
    result$bbox_overlaps == 0L, result$leader_crossings == 0L)
  results[[name]] <- result
}
print(do.call(rbind, results), row.names = FALSE)
