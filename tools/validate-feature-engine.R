# Explicit visual acceptance for the v0.13.1 circular Feature Engine.
# Run from the package root with:
#   Rscript tools/validate-feature-engine.R [output-directory]
#
# The 13 checked-in reference FASTA/.dna/PNG trios define the visual benchmark.
# Synthetic class fixtures below remain deliberately generic and must never be
# treated as exact records of a named vector.

devtools::load_all(quiet = TRUE)

crossing_probe <- data.frame(x1 = 0, y1 = 0, x2 = 2, y2 = 0)
occupied_probe <- matrix(c(.5, -1, .5, 1, 1.5, -1, 1.5, 1),
  nrow = 2L, byrow = TRUE)
stopifnot(ggchord:::ggchord_external_route_crossing_count(
  crossing_probe, occupied_probe
) == 2L)

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

requested_map_arg <- commandArgs(trailingOnly = TRUE)[2L]
requested_names <- if (!is.na(requested_map_arg) && nzchar(requested_map_arg)) {
  requested_maps <- trimws(strsplit(requested_map_arg, ",", fixed = TRUE)[[1L]])
  unique(c(requested_maps, make.names(requested_maps)))
} else names(reference_maps)
selected_reference_maps <- reference_maps[
  names(reference_maps) %in% requested_names |
    unname(reference_maps) %in% requested_names
]

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
reference_restriction_sites <- function(sequence, name) {
  sites <- find_restriction_sites(sequence)
  if (identical(name, "pETDuet-1")) {
    # The saved .dna profile is "None", but the checked-in PNG was exported
    # with a dense unique-6+ display. The visual fixture therefore owns this
    # explicit override; it must not be presented as recovered file metadata.
    return(filter_restriction_sites(sites, set = "unique_6plus",
      parent_set = "commercial_nonredundant"))
  }
  filter_restriction_sites(sites, set = "reference",
    parent_set = "commercial_nonredundant")
}
real_fixtures <- lapply(unname(selected_reference_maps), function(name) {
  sequence <- read_reference_fasta(name)
  sites <- reference_restriction_sites(sequence, name)
  primers <- find_primer_bindings(sequence, set = "reference")
  list(
    sequence = sequence, features = find_common_features(sequence),
    sites = sites, primers = primers
  )
})
names(real_fixtures) <- names(selected_reference_maps)
if (identical(names(selected_reference_maps), names(reference_maps))) {
  stopifnot(
    nrow(real_fixtures[[make.names("pSpCas9(BB)-2A-GFP (PX458)")]]$sites) == 4L,
    nrow(real_fixtures[[make.names("pTRIPZ")]]$sites) == 38L,
    nrow(real_fixtures[[make.names("pETDuet-1")]]$sites) == 58L,
    sum(vapply(real_fixtures, function(x) nrow(x$primers), integer(1L))) == 7L
  )
}
if ("pSB1C3" %in% names(real_fixtures)) {
  psb_sites <- real_fixtures[["pSB1C3"]]$sites
  stopifnot(any(psb_sites$enzyme == "PflMI" &
    psb_sites$display_warning == "methylation_blocked"))
}
if (make.names("pSpCas9(BB)-2A-GFP (PX458)") %in% names(real_fixtures)) {
  px_features <- real_fixtures[[make.names(
    "pSpCas9(BB)-2A-GFP (PX458)")]]$features
  stopifnot(px_features$feature_shape[px_features$anno == "hybrid intron"] ==
    "capped_line")
}
if ("pTRE.Tight.BI" %in% names(real_fixtures)) {
  tre_features <- real_fixtures[["pTRE.Tight.BI"]]$features
  stopifnot(tre_features$strand[tre_features$anno ==
    "bidirectional TRE promoter"] == "+/-")
}
bidirectional_promoter <- ggchord:::ggchord_feature_geometry(
  "promoter_arrow", 0, .08, .9, .07, 1, ref = list(),
  draw_head = TRUE, draw_start_head = TRUE, bidirectional = TRUE
)[[1L]]
stopifnot(sum(abs(bidirectional_promoter$radius - .9) < 1e-10) == 2L)
axis_lengths <- c(2070, 2686, 4361, 5369, 9288, 13320)
axis_steps <- vapply(axis_lengths, function(x) {
  diff(ggchord:::breakPointsFunc(x))[1L]
}, numeric(1L))
stopifnot(identical(axis_steps, c(250, 500, 500, 1000, 1000, 2000)))
reference_fixture_names <- names(real_fixtures)
fixtures <- c(real_fixtures, fixtures)

fixture_names <- names(fixtures)
if (!is.na(commandArgs(trailingOnly = TRUE)[2L]) &&
    nzchar(commandArgs(trailingOnly = TRUE)[2L])) {
  fixture_names <- fixture_names[fixture_names %in% requested_names]
  if (!length(fixture_names)) stop("No requested validation map was found")
}

for (name in fixture_names) {
  message("Validating ", name, "...")
  fixture <- fixtures[[name]]
  tracks <- position_feature_stack(
    spacing = .085, base_position = position_plasmid()
  )
  primer_layers <- if (!is.null(fixture$primers) && nrow(fixture$primers)) {
    list(
      geom_primer(data = fixture$primers),
      geom_primer_label_repel(data = fixture$primers, max_overlaps = 0)
    )
  } else NULL
  plot <- ggchord(fixture$sequence, validate = "none") +
    geom_seq(seq_style = "double") +
    geom_feature_plasmid(data = fixture$features, position = tracks) +
    geom_feature_label_repel(data = fixture$features, position = tracks,
      external = TRUE, max_overlaps = 0) +
    (if (!is.null(fixture$sites)) {
      geom_restriction_site(data = fixture$sites)
    } else NULL) +
    primer_layers +
    geom_seq_center_label() + coord_circular(rotation = 90) +
    theme_ggchord_plasmid()
  # Every reference must remain exportable on a constrained square device;
  # automatic preview below then rebuilds from a clean cache at its intended
  # content-derived size.
  close_device <- ggchord:::ggchord_measurement_device(width = 6, height = 6)
  layout <- tryCatch(
    export_ggchord_layout(
      plot, include = c("feature", "labels", "restriction")
    ),
    finally = close_device()
  )
  message("Exported geometry for ", name)
  stopifnot(all(is.finite(layout$feature$x)), all(is.finite(layout$feature$y)))
  if (identical(unname(reference_maps[name]), "pTRIPZ")) stopifnot(
    sum(layout$labels$.component == "text" &
      layout$labels$anno == "tet operator" &
      layout$labels$feature_label_mode == "external", na.rm = TRUE) == 6L
  )
  if (identical(unname(reference_maps[name]), "pTRE-Tight-BI")) stopifnot(
    sum(layout$labels$.component == "text" &
      layout$labels$anno == "tet operator" &
      layout$labels$feature_label_mode == "external", na.rm = TRUE) == 0L
  )
  inner_text <- layout$labels[
    layout$labels$.component == "text" &
      layout$labels$feature_label_mode %in% c("inside", "adjacent") &
      is.finite(layout$labels$feature_track) &
      is.finite(layout$labels$label_track), , drop = FALSE
  ]
  for (sid in unique(inner_text$accver)) {
    rows <- which(inner_text$accver == sid)
    feature_tracks <- sort(unique(inner_text$feature_track[rows]))
    for (row in rows) {
      next_track <- feature_tracks[
        feature_tracks > inner_text$feature_track[row]
      ][1L]
      if (is.finite(next_track)) stopifnot(
        inner_text$label_track[row] < next_track
      )
    }
  }
  stopifnot(all(c("annotation_class", "dominant_band", "spill_reason",
    "leader_crossing_count") %in% names(layout$annotation_registry)))
  visible_boxes <- layout$annotation_registry[
    is.finite(layout$annotation_registry$bbox_xmin) &
      is.finite(layout$annotation_registry$bbox_xmax) &
      is.finite(layout$annotation_registry$bbox_ymin) &
      is.finite(layout$annotation_registry$bbox_ymax), , drop = FALSE
  ]
  if (nrow(visible_boxes)) stopifnot(
    all(visible_boxes$bbox_xmin >= layout$metadata$xlim[1L] - 1e-5),
    all(visible_boxes$bbox_xmax <= layout$metadata$xlim[2L] + 1e-5),
    all(visible_boxes$bbox_ymin >= layout$metadata$ylim[1L] - 1e-5),
    all(visible_boxes$bbox_ymax <= layout$metadata$ylim[2L] + 1e-5)
  )
  measured_crossings <- layout$annotation_registry$leader_crossing_count
  positive_crossings <- measured_crossings[is.finite(measured_crossings) &
    measured_crossings > 0L]
  leader_error <- NULL
  if (length(positive_crossings)) {
    bad <- layout$annotation_registry[
      is.finite(measured_crossings) & measured_crossings > 0L, , drop = FALSE
    ]
    bad_labels <- layout$restriction[
      layout$restriction$restriction_component == "label" &
        layout$restriction$source_row %in% bad$source_row,
      c("source_row", "anchor_position", "label"), drop = FALSE
    ]
    leader_error <- paste0(
      as.character(fixture$sequence$label[1L]),
      " has ", length(positive_crossings),
      " registry rows on crossing leaders (max ", max(positive_crossings),
      "): ", paste(paste0(bad_labels$source_row, "@",
        bad_labels$anchor_position, " ", bad_labels$label), collapse = "; ")
    )
  }
  restriction_registry <- layout$annotation_registry[
    !is.na(layout$annotation_registry$annotation_class) &
      layout$annotation_registry$annotation_class == "restriction", ,
    drop = FALSE
  ]
  if (nrow(restriction_registry)) {
    used_bands <- sort(unique(restriction_registry$band))
    stopifnot(identical(used_bands, seq_len(max(used_bands))))
  }
  restriction_labels <- layout$restriction[
    layout$restriction$restriction_component == "label", , drop = FALSE
  ]
  sequence_label <- as.character(fixture$sequence$label[1L])
  if (sequence_label == "pBR322") {
    bottom <- restriction_labels$anchor_position >= 2060 &
      restriction_labels$anchor_position <= 2360
    stopifnot(
      sum(bottom) >= 8L,
      !any(restriction_labels$label_layout[bottom] == "perimeter_rail"),
      any(restriction_labels$outer_track[bottom] > 1L)
    )
  }
  if (sequence_label == "pUC19") {
    stopifnot(
      !any(restriction_labels$label_layout == "perimeter_rail"),
      mean(restriction_labels$outer_track == 1L) >= .80
    )
  }
  if (sequence_label == "pBluescript II SK(+)") {
    right <- restriction_labels$anchor_position >= 650 &
      restriction_labels$anchor_position <= 760
    upper_left <- restriction_labels$anchor_position >= 2520 &
      restriction_labels$anchor_position <= 2650
    stopifnot(
      sum(right) >= 12L,
      !any(restriction_labels$label_layout[right] == "perimeter_rail"),
      !any(restriction_labels$label_layout[upper_left] == "perimeter_rail")
    )
  }
  if (sequence_label == "pETDuet-1") {
    dense <- restriction_labels$anchor_position <= 451
    stopifnot(
      sum(dense) >= 10L,
      !any(restriction_labels$label_layout == "perimeter_rail"),
      mean(restriction_labels$outer_track[dense] == 1L) >= .80
    )
  }
  text <- layout$labels[layout$labels$.component == "text", , drop = FALSE]
  if (identical(as.character(fixture$sequence$label[1L]),
      "pBluescript II SK(+)")) {
    raw_layout <- get_chord_layout(plot)
    major_ticks <- raw_layout$axis_ticks[raw_layout$axis_ticks$is_major, ]
    labelled_ticks <- major_ticks[!is.na(major_ticks$label), ]
    stopifnot(
      all(abs(sqrt(major_ticks$x0^2 + major_ticks$y0^2) - .9875) < 1e-6),
      all(labelled_ticks$label_along_axis),
      all(labelled_ticks$label_hjust == 0),
      all(labelled_ticks$label_vjust == .5)
    )
    site_ticks <- layout$restriction[
      layout$restriction$restriction_component == "tick", , drop = FALSE
    ]
    site_roots <- site_ticks[!duplicated(site_ticks$group), , drop = FALSE]
    stopifnot(all(abs(sqrt(site_roots$x^2 + site_roots$y^2) - 1.0125) < 3e-5))
    feature_lanes <- unique(layout$feature[
      layout$feature$.component == "polygon", c("anno", "lane")
    ])
    lane_of <- stats::setNames(feature_lanes$lane, feature_lanes$anno)
    stopifnot(
      all(lane_of[c(
        "lacZα", "lac operator", "lac promoter", "ori", "AmpR",
        "AmpR promoter"
      )] == 0L),
      all(lane_of[c(
        "f1 ori", "M13 fwd", "T7 promoter", "MCS", "T3 promoter", "M13 rev"
      )] == 1L),
      all(lane_of[c("KS primer", "SK primer")] == 2L)
    )
    leader_features <- unique(layout$labels$anno[
      layout$labels$.component == "segment"
    ])
    expected_leaders <- text$anno[
      text$label_track > text$feature_track + 1L
    ]
    # A visibly wide immediate gutter may itself need a short connector;
    # long cross-track moves must have one, but adjacent-track leaders are
    # permitted when their real radial gap warrants it.
    stopifnot(all(expected_leaders %in% leader_features))
    stopifnot(all(leader_features %in% text$anno[
      text$feature_label_mode %in% c("adjacent", "external")
    ]))
    boundary_features <- table(layout$feature$anno[
      layout$feature$.component == "boundary"
    ])
    stopifnot(
      unname(boundary_features["AmpR"]) == 8L,
      unname(boundary_features["lac promoter"]) == 16L
    )
    polygon <- layout$feature[
      layout$feature$.component == "polygon", , drop = FALSE
    ]
    polygon$.radius <- sqrt(polygon$x^2 + polygon$y^2)
    band_bounds <- lapply(split(polygon, polygon$lane), function(x) {
      range(x$.radius)
    })
    stopifnot(
      band_bounds[["0"]][1L] > band_bounds[["1"]][2L],
      band_bounds[["1"]][1L] > band_bounds[["2"]][2L]
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
    arc_text <- text[
      text$feature_label_mode != "external", , drop = FALSE
    ]
    stopifnot(
      all(arc_text$.draw_as_arc),
      nrow(glyphs) == sum(nchar(arc_text$label, type = "chars"))
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
  external_text <- feature_external$label
  if (nrow(feature_external) && "annotation_class" %in% names(feature_external) &&
      "feature_label" %in% names(feature_external)) {
    primer_text <- !is.na(feature_external$annotation_class) &
      feature_external$annotation_class == "primer" &
      !is.na(feature_external$feature_label)
    external_text[primer_text] <- feature_external$feature_label[primer_text]
  }
  external <- if (nrow(feature_external)) data.frame(
      text_x = feature_external$x, text_y = feature_external$y,
      text = external_text, size = feature_external$size,
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
    if (any(collisions)) {
      leader_error <- paste0(
        leader_error %||% "", " External feature/primer labels overlap in ",
        as.character(fixture$sequence$label[1L]), ": ",
        paste(which(collisions), collapse = ", ")
      )
    }
  }
  if (name %in% reference_fixture_names) {
    plot$ggchord$ref$layout <- NULL
    preview <- view_ggchord(plot, viewer = "none")
    output_file <- file.path(output_dir, paste0(name, ".png"))
    stopifnot(file.copy(preview, output_file, overwrite = TRUE))
    if (nrow(fixture$primers)) {
      primer_registry <- layout$annotation_registry[
        !is.na(layout$annotation_registry$annotation_class) &
          layout$annotation_registry$annotation_class == "primer", ,
        drop = FALSE
      ]
      stopifnot(nrow(primer_registry) == nrow(fixture$primers))
      primer_text <- layout$labels[
        layout$labels$.component == "text" &
          !is.na(layout$labels$annotation_class) &
          layout$labels$annotation_class == "primer", , drop = FALSE
      ]
      stopifnot(
        all(fixture$primers$name %in% primer_text$primer_name),
        all(primer_text$primer_show_location),
        all(grepl("\\([0-9]+ \\.\\. [0-9]+\\)", primer_text$feature_label))
      )
    }
  } else {
    ggplot2::ggsave(file.path(output_dir, paste0(name, ".png")), plot,
      width = 6, height = 6, dpi = 144)
  }
  if (!is.null(leader_error)) stop(leader_error)
}

if (is.na(requested_map_arg) || !nzchar(requested_map_arg)) {
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
}

writeLines(c(
  "Reference maps inspected from examples/plasmid:",
  unname(selected_reference_maps),
  paste0("Rendered fixtures: ", output_dir)
))
