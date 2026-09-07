feature_test_data <- function() {
  data.frame(
    accver = c("A", "A", "B", "B"),
    start = c(100, 300, 100, 300), end = c(220, 420, 220, 420),
    strand = rep(c("+", "-"), 2), anno = letters[1:4],
    stringsAsFactors = FALSE
  )
}

feature_position_export <- function(position = "identity", geom = "gene",
                                    data = feature_test_data()) {
  seq <- data.frame(accver = c("A", "B"), length = 1000)
  layer <- if (geom == "gene") {
    geom_gene(data = data, position = position)
  } else {
    geom_feature(data = data, position = position, feature_shape = "arrow")
  }
  category <- if (geom == "gene") "gene" else "feature"
  export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() + layer,
    include = category
  )[[category]]
}

test_that("gene and feature identity really use the sequence centerline", {
  gene <- feature_position_export()
  feature <- feature_position_export(geom = "feature")
  expect_equal(unique(gene$position_name), "identity")
  expect_equal(unique(gene$normal_offset), 0)
  expect_equal(unique(feature$normal_offset), 0)
})

test_that("identity Position accepts ggplot2 string, constructor, and prototype", {
  string <- feature_position_export("identity")
  constructor <- feature_position_export(position_identity())
  prototype <- feature_position_export(PositionIdentity)
  columns <- c("x", "y", "source_row", "normal_offset")
  expect_equal(string[columns], constructor[columns])
  expect_equal(string[columns], prototype[columns])
  expect_true("position_identity" %in% getNamespaceExports("ggchord"))
  expect_true("PositionIdentity" %in% getNamespaceExports("ggchord"))
})

test_that("strand Position implements the signed outward contract", {
  standard <- feature_position_export(position_strand(.1))
  offsets <- vapply(split(standard$normal_offset, standard$source_row), unique, numeric(1))
  expect_equal(unname(offsets), c(.1, -.1, .1, -.1))

  reversed <- feature_position_export(position_strand(-.1))
  offsets <- vapply(split(reversed$normal_offset, reversed$source_row), unique, numeric(1))
  expect_equal(unname(offsets), c(-.1, .1, -.1, .1))

  named <- feature_position_export(position_strand(c("+" = -.08, "-" = .05)))
  offsets <- vapply(split(named$normal_offset, named$source_row), unique, numeric(1))
  expect_equal(unname(offsets), c(-.08, .05, -.08, .05))
})

test_that("flexible sequence and strand Position specifications survive", {
  by_sequence <- feature_position_export(position_strand(c(A = .08, B = .14)))
  offsets <- vapply(split(by_sequence$normal_offset, by_sequence$source_row), unique, numeric(1))
  expect_equal(unname(offsets), c(.08, -.08, .14, -.14))

  combined <- feature_position_export(position_strand(list(
    A = c("+" = .07, "-" = -.04),
    B = c("+" = -.12, "-" = .03)
  )))
  offsets <- vapply(split(combined$normal_offset, combined$source_row), unique, numeric(1))
  expect_equal(unname(offsets), c(.07, -.04, -.12, .03))
})

test_that("plasmid Position shares one band without changing strand", {
  plasmid <- feature_position_export("plasmid")
  expect_equal(unique(plasmid$normal_offset), -.1)
  expect_equal(sort(unique(plasmid$strand)), c("-", "+"))

  x <- -.15
  shortcut <- feature_position_export(position_plasmid(x))
  explicit <- feature_position_export(position_strand(c("+" = x, "-" = x)))
  expect_equal(shortcut$x, explicit$x)
  expect_equal(shortcut$y, explicit$y)
  expect_equal(shortcut$source_row, explicit$source_row)

  zero <- feature_position_export(position_plasmid(0))
  identity <- feature_position_export("identity")
  expect_equal(zero$x, identity$x)
  expect_equal(zero$y, identity$y)
})

test_that("legacy gene_offset translates exactly and conflicts clearly", {
  seq <- data.frame(accver = c("A", "B"), length = 1000)
  data <- feature_test_data()
  legacy_plot <- expect_warning(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_gene(data = data, gene_offset = .1),
    "deprecated"
  )
  legacy <- export_ggchord_layout(legacy_plot, include = "gene")$gene
  modern <- feature_position_export("strand")
  expect_equal(legacy$x, modern$x)
  expect_equal(legacy$y, modern$y)

  flexible_plot <- expect_warning(
    ggchord(seq, validate = "none") + geom_seq() + geom_gene(
      data = data,
      gene_offset = list(
        A = c("+" = .07, "-" = .04),
        B = c("+" = .12, "-" = .03)
      )
    ),
    "deprecated"
  )
  flexible <- export_ggchord_layout(flexible_plot, include = "gene")$gene
  flexible_offsets <- vapply(
    split(flexible$normal_offset, flexible$source_row), unique, numeric(1)
  )
  expect_equal(unname(flexible_offsets), c(.07, -.04, .12, -.03))

  expect_error(
    geom_gene(gene_offset = .1, position = "plasmid"), "cannot be combined"
  )
  expect_error(
    geom_feature(feature_offset = .1, position = "strand"), "cannot be combined"
  )
})

test_that("feature stack composes after its base Position", {
  data <- data.frame(
    accver = "A", start = c(100, 150, 500), end = c(250, 230, 600),
    strand = c("+", "-", "+"), anno = letters[1:3]
  )
  identity <- feature_position_export(
    position_feature_stack(base_position = "identity"), data = data
  )
  info <- unique(identity[c("source_row", "base_offset", "lane", "normal_offset")])
  expect_equal(info$base_offset, rep(0, nrow(info)))
  expect_equal(info$lane[order(info$source_row)], c(0L, 1L, 0L))

  strand <- feature_position_export(
    position_feature_stack(base_position = "strand"), data = data
  )
  strand_info <- unique(strand[c("source_row", "base_offset", "lane")])
  expect_equal(strand_info$base_offset[order(strand_info$source_row)], c(.1, -.1, .1))
  expect_equal(strand_info$lane, rep(0L, nrow(strand_info)))

  plasmid <- feature_position_export(
    position_feature_stack(base_position = "plasmid"), data = data
  )
  plasmid_info <- unique(plasmid[c("source_row", "base_offset", "lane")])
  expect_equal(plasmid_info$base_offset, rep(-.1, nrow(plasmid_info)))
  expect_equal(plasmid_info$lane[order(plasmid_info$source_row)], c(0L, 1L, 0L))

  legacy <- feature_position_export(position_feature_stack(), data = data)
  expect_equal(unique(legacy$position_name), "feature_stack_legacy")

  overlapping <- data.frame(
    accver = "A", start = c(100, 150, 100, 150),
    end = c(250, 230, 250, 230), strand = c("+", "+", "-", "-"),
    anno = letters[1:4]
  )
  legacy_lanes <- feature_position_export(
    position_feature_stack(), data = overlapping
  )
  legacy_info <- unique(legacy_lanes[
    c("source_row", "lane", "normal_offset")
  ])
  legacy_info <- legacy_info[order(legacy_info$source_row), ]
  expect_equal(legacy_info$lane, c(0L, 1L, 0L, 1L))
  expect_equal(legacy_info$normal_offset, c(.1, .18, -.1, -.18))
  expect_error(
    position_feature_stack(side = "inside", base_position = "identity"),
    "cannot be combined"
  )
})

test_that("signed Position metadata is consistent across coordinates", {
  seq <- data.frame(accver = "A", length = 1000)
  data <- feature_test_data()[1:2, ]
  chord <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_gene(data = data, position = "strand") + coord_chord(rotation = 0),
    include = "gene"
  )$gene
  circle <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_gene(data = data, position = "strand") + coord_circular(),
    include = "gene"
  )$gene
  expect_equal(unique(chord[c("source_row", "normal_offset")]),
    unique(circle[c("source_row", "normal_offset")]))
  expect_true(all(c(
    "position_name", "base_offset", "lane", "lane_offset", "normal_offset"
  ) %in% names(circle)))
})
