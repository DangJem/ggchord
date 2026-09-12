# Explicit visual acceptance for the v0.13.0 circular Feature Engine.
# Run from the package root with:
#   Rscript tools/validate-feature-engine.R [output-directory]
#
# Reference maps are inspected manually; their copyrighted .dna annotations
# are not redistributed here. Synthetic class fixtures below are deliberately
# named as classes and must never be treated as exact records of those vectors.

devtools::load_all(quiet = TRUE)

reference_maps <- c(
  "pBluescript II SK(+)" = "https://www.snapgene.com/plasmids/basic_cloning_vectors/pBluescript_II_SK%28%2B%29",
  "pUC19 / pUC19c" = "https://www.snapgene.com/plasmids/basic_cloning_vectors/pUC19",
  "pET-28a(+)" = "https://www.snapgene.com/plasmids/pet_and_duet_vectors_%28novagen%29/pET-28a%28%2B%29",
  "pLKO.1" = "https://www.snapgene.com/plasmids/viral_expression_and_packaging_vectors/pLKO.1",
  "pLEX-MCS" = "https://www.snapgene.com/plasmids/viral_expression_and_packaging_vectors/pLEX-MCS",
  "pPICZ(alpha) A" = "https://www.snapgene.com/plasmids/yeast_plasmids/pPICZ%28alpha%29_A",
  "pcDNA3.1 CT-GFP" = "https://www.snapgene.com/plasmids/mammalian_expression_vectors/pcDNA3.1_CT-GFP",
  "pEGFP-N1" = "https://www.snapgene.com/plasmids/fluorescent_protein_genes_and_plasmids/pEGFP-N1"
)

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

data(plasmid_example_pBluescript_II_SK_plus)
data(plasmid_example_pUC19c)
real_fixtures <- list(
  pBluescript_II_SK_plus = list(
    sequence = plasmid_example_pBluescript_II_SK_plus,
    features = find_common_features(plasmid_example_pBluescript_II_SK_plus)
  ),
  pUC19c = list(
    sequence = plasmid_example_pUC19c,
    features = find_common_features(plasmid_example_pUC19c)
  )
)
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
    geom_seq_center_label() + coord_circular(rotation = 90) +
    theme_ggchord_plasmid()
  layout <- export_ggchord_layout(plot, include = c("feature", "labels"))
  stopifnot(all(is.finite(layout$feature$x)), all(is.finite(layout$feature$y)))
  text <- layout$labels[layout$labels$.component == "text", , drop = FALSE]
  external <- text[text$feature_label_mode == "external", , drop = FALSE]
  if (nrow(external) > 1L) {
    boxes <- ggchord:::ggchord_text_boxes(
      external, units_per_inch = .35, box_padding = .015
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
  "Reference maps inspected manually (not redistributed):",
  paste(names(reference_maps), unname(reference_maps), sep = "\t"),
  paste0("Rendered fixtures: ", output_dir)
))
