# Explicit visual acceptance for the v0.13.0 circular Feature Engine.
# Run from the package root with:
#   Rscript tools/validate-feature-engine.R [output-directory]
#
# The 13 checked-in reference FASTA/.dna/PNG trios define the visual benchmark.
# Synthetic class fixtures below remain deliberately generic and must never be
# treated as exact records of a named vector.

devtools::load_all(quiet = TRUE)

reference_maps <- c(
  "pBR322", "pUC19", "pBluescript II SK(+)", "pSB1C3", "pET-28a(+)",
  "pETDuet-1", "pcDNA3.1(+)", "pTRE-Tight-BI",
  "pSpCas9(BB)-2A-GFP (PX458)", "pDONR221", "pCAMBIA1300",
  "pEarleyGate 201", "pTRIPZ"
)
names(reference_maps) <- make.names(reference_maps)
reference_paths <- lapply(c(".fna", ".dna", ".png"), function(extension) {
  file.path("examples", "plasmid", paste0(unname(reference_maps), extension))
})
stopifnot(all(vapply(reference_paths, function(path) all(file.exists(path)),
  logical(1L))))

output_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(output_dir) || !nzchar(output_dir)) {
  output_dir <- file.path(tempdir(), "ggchord-feature-engine")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

class_fixture <- function(name, length, starts, ends, directionality, labels,
                          priority = 0) {
  list(
    sequence = data.frame(accver = name, label = name, length = length),
    features = data.frame(
      accver = name, start = starts, end = ends,
      directionality = directionality, anno = labels,
      display_priority = priority, stringsAsFactors = FALSE
    )
  )
}

fixtures <- list(
  bacterial_expression_class = class_fixture(
    "bacterial-expression-class", 5400,
    c(80, 120, 155, 190, 225, 900, 2100, 3500),
    c(560, 205, 190, 250, 310, 1800, 2900, 4400),
    c("forward", "forward", "nondirectional", "forward", "forward",
      "reverse", "forward", "nondirectional"),
    c("long control", "short tag", "point-like motif", "short operator",
      "short binding region", "long coding feature", "marker feature",
      "origin feature")
  ),
  lentiviral_class = class_fixture(
    "lentiviral-class", 10500,
    c(50, 180, 700, 1150, 1900, 3100, 3300, 6100, 7300, 9000),
    c(650, 900, 1700, 2600, 3500, 3550, 4700, 7200, 8400, 10200),
    c("forward", "forward", "nondirectional", "forward", "forward",
      "bidirectional", "forward", "reverse", "forward", "forward"),
    paste("viral feature", seq_len(10))
  ),
  mammalian_expression_class = class_fixture(
    "mammalian-expression-class", 6100,
    c(40, 220, 410, 570, 650, 1450, 1700, 2600, 3900, 4700),
    c(500, 620, 650, 720, 1550, 1680, 2500, 3300, 4650, 5750),
    c("nondirectional", "forward", "forward", "nondirectional", "forward",
      "nondirectional", "forward", "forward", "reverse", "nondirectional"),
    paste("mammalian feature", seq_len(10))
  ),
  yeast_expression_class = class_fixture(
    "yeast-expression-class", 3600,
    c(2, 900, 940, 1010, 1060, 1180, 1420, 1700, 2100, 2600),
    c(940, 1260, 1210, 1040, 1080, 1360, 1660, 2090, 2540, 3250),
    c("forward", "forward", "forward", "forward", "forward", "nondirectional",
      "nondirectional", "forward", "forward", "nondirectional"),
    paste("yeast feature", seq_len(10))
  )
)

stress_segments <- data.frame(
  segment_index = 1:4,
  segment_type = c("standard", "standard", "gap", "standard"),
  start = c(600, 681, 741, 781), end = c(680, 740, 780, 920),
  color = c("#D95F02", "#7570B3", NA, "#1B9E77"),
  segment_name = c("one", "two", NA, "three")
)
stress <- class_fixture(
  "synthetic-stress", 2000,
  c(100, 100, 160, 240, 300, 300, 330, 600, 1050, 1200, 1400),
  c(520, 520, 360, 480, 620, 620, 370, 920, 1050, 1206, 1850),
  c("forward", "reverse", "nondirectional", "forward", "forward", "reverse",
    "forward", "forward", "nondirectional", "forward", "bidirectional"),
  c("complete overlap A", "complete overlap B", "containment", "chain",
    "same length A", "same length B", "prioritized short", "segmented",
    "point", "very short", "a deliberately very long external label"),
  priority = c(rep(0, 6), 10, rep(0, 4))
)
stress$features$segments <- vector("list", nrow(stress$features))
stress$features$segments[[8]] <- stress_segments
fixtures$synthetic_stress <- stress

read_reference_fasta <- function(name) {
  path <- file.path("examples", "plasmid", paste0(name, ".fna"))
  lines <- readLines(path, warn = FALSE)
  sequence <- paste0(lines[!grepl("^>", lines)], collapse = "")
  data.frame(
    accver = make.names(name), label = name, length = nchar(sequence),
    sequence = sequence, stringsAsFactors = FALSE
  )
}
real_fixtures <- lapply(unname(reference_maps), function(name) {
  sequence <- read_reference_fasta(name)
  sites <- find_restriction_sites(sequence)
  sites <- filter_restriction_sites(
    sites, set = "unique_6plus", parent_set = "commercial_nonredundant"
  )
  list(
    sequence = sequence, features = find_common_features(sequence),
    sites = sites
  )
})
names(real_fixtures) <- names(reference_maps)
fixtures <- c(real_fixtures, fixtures)

for (name in names(fixtures)) {
  fixture <- fixtures[[name]]
  tracks <- position_feature_stack(
    spacing = .085, base_position = position_plasmid()
  )
  plot <- ggchord(fixture$sequence, validate = "none") +
    geom_seq(seq_style = "double") +
    geom_feature_plasmid(data = fixture$features, position = tracks) +
    geom_feature_label_repel(data = fixture$features, position = tracks,
      external = TRUE, max_overlaps = 0) +
    (if (!is.null(fixture$sites)) {
      geom_restriction_site(data = fixture$sites)
    } else NULL) +
    geom_seq_center_label() + coord_circular(rotation = 90) +
    theme_ggchord_plasmid()
  layout <- export_ggchord_layout(
    plot, include = c("feature", "labels", "restriction")
  )
  stopifnot(all(is.finite(layout$feature$x)), all(is.finite(layout$feature$y)))
  text <- layout$labels[layout$labels$.component == "text", , drop = FALSE]
  if (identical(as.character(fixture$sequence$label[1L]),
      "pBluescript II SK(+)")) {
    feature_lanes <- unique(layout$feature[
      layout$feature$.component == "polygon", c("anno", "lane")
    ])
    lane_of <- stats::setNames(feature_lanes$lane, feature_lanes$anno)
    stopifnot(
      all(lane_of[c(
        "M13 fwd", "T7 promoter", "MCS", "T3 promoter", "M13 rev",
        "lac operator", "lac promoter"
      )] == 1L),
      all(lane_of[c("KS primer", "SK primer")] == 2L),
      !any(text$feature_label_mode == "external")
    )
    polygon <- layout$feature[
      layout$feature$.component == "polygon", , drop = FALSE
    ]
    polygon$.radius <- sqrt(polygon$x^2 + polygon$y^2)
    band_bounds <- lapply(split(polygon, polygon$lane), function(x) {
      range(x$.radius)
    })
    stopifnot(
      band_bounds[["0"]][1L] - band_bounds[["1"]][2L] > .10,
      band_bounds[["1"]][1L] - band_bounds[["2"]][2L] > .10
    )
    adjacent <- text[text$feature_label_mode == "adjacent", , drop = FALSE]
    radius <- sqrt(adjacent$x^2 + adjacent$y^2)
    stopifnot(all(tapply(radius, adjacent$label_track,
      function(value) diff(range(value))) < 1e-8))
    tangent <- atan2(adjacent$y, adjacent$x) * 180 / pi + 90
    error <- abs((adjacent$angle - tangent + 90) %% 180 - 90)
    stopifnot(all(error < 1e-6))
    glyphs <- layout$labels[
      layout$labels$.component == "arc_text", , drop = FALSE
    ]
    stopifnot(
      all(text$.draw_as_arc),
      nrow(glyphs) == sum(nchar(text$label, type = "chars"))
    )
    glyph_radius <- sqrt(glyphs$x^2 + glyphs$y^2)
    glyph_tangent <- atan2(glyphs$y, glyphs$x) * 180 / pi + 90
    stopifnot(
      all(tapply(glyph_radius, glyphs$.arc_parent_group,
        function(value) diff(range(value))) < 1e-8),
      all(abs((glyphs$angle - glyph_tangent + 90) %% 180 - 90) < 1e-6)
    )
  }
  feature_external <- text[
    text$feature_label_mode == "external", , drop = FALSE
  ]
  # Restriction labels own their circular contour/fan and must not be folded
  # into the feature-callout collision rail. Validate feature callouts only;
  # restriction geometry has its own ordered-fan checks.
  external <- if (nrow(feature_external)) data.frame(
      text_x = feature_external$x, text_y = feature_external$y,
      text = feature_external$label, size = feature_external$size,
      hjust = feature_external$hjust, vjust = feature_external$vjust,
      text_angle = 0
    ) else data.frame()
  if (nrow(external) > 1L) {
    boxes <- ggchord:::ggchord_text_boxes(
      external, units_per_inch = .35, box_padding = .015,
      x_col = "text_x", y_col = "text_y", text_col = "text",
      angle_col = "text_angle"
    )
    pairs <- utils::combn(seq_len(nrow(boxes)), 2L)
    collisions <- apply(pairs, 2L, function(pair) {
      ggchord:::ggchord_oriented_box_overlaps(
        boxes[pair[1L], , drop = FALSE],
        boxes[pair[2L], , drop = FALSE]
      )
    })
    stopifnot(!any(collisions))
  }
  ggplot2::ggsave(file.path(output_dir, paste0(name, ".png")), plot,
    width = 6, height = 6, dpi = 144)
}

stress_layout <- export_ggchord_layout(
  ggchord(stress$sequence, validate = "none") + geom_seq() +
    geom_feature_plasmid(data = stress$features) + coord_circular(),
  include = "feature"
)$feature
stress_lanes <- unique(stress_layout[, c("anno", "lane")])
stopifnot(
  stress_lanes$lane[stress_lanes$anno == "prioritized short"] == 0L,
  unique(stress_layout$.component[stress_layout$anno == "point"]) == "point",
  length(unique(stress_layout$lane[stress_layout$anno == "segmented"])) == 1L
)

writeLines(c(
  "Reference maps inspected from examples/plasmid:",
  unname(reference_maps),
  paste0("Rendered fixtures: ", output_dir)
))
