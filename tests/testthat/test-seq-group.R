test_that("geom_seq remains backward compatible as a single layer", {
  layers <- geom_seq()
  expect_length(layers, 1)
  expect_identical(layers[[1]]$ggchord_type, "seq")
})

test_that("seq_group from a seq_data column draws group labels and adds a group gap", {
  data(seq_data_example)
  seq_data_example$seq_group <- c("even", "even", "odd", "odd")

  p <- ggchord(seq_data_example, validate = "none") +
    geom_seq(seq_group_gap = 0.08)

  # geom_seq() itself still returns a single layer
  expect_length(geom_seq(), 1)

  built <- ggplot_build(p)
  expect_s3_class(built, "ggplot_built")

  layout <- get_chord_layout()
  expect_equal(sort(unique(layout$seq_groups)), c("even", "odd"))
  expect_equal(nrow(layout$group_labels), 2)
  expect_equal(sort(layout$group_labels$label), c("even", "odd"))
  expect_true(all(layout$group_labels$text_angle == 0))

  # Group gap widens the boundary between the first two differently-grouped
  # sequences relative to an ungrouped plot.
  p0 <- ggchord(seq_data_example[, c("seq_id", "length")], validate = "none") +
    geom_seq()
  ggplot_build(p0)
  layout0 <- get_chord_layout()

  seqs <- layout$seqs
  arc_start <- function(layout, id) {
    idx <- which(vapply(layout$seq_arcs, function(arc) identical(unique(arc$seq_id)[1], id), logical(1)))
    arc <- layout$seq_arcs[[idx]]
    atan2(arc$y[1], arc$x[1])
  }
  expect_equal(seqs, layout0$seqs)

  # First boundary is odd -> even in the grouped data; the grouped start of
  # the second sequence must be larger than the ungrouped start.
  expect_gt(arc_start(layout, seqs[3]), arc_start(layout0, seqs[3]))
})

test_that("seq_group accepts named vectors and per-group labels/colors", {
  data(seq_data_example)
  grp <- c("MT108731.1" = "host", "MT118296.1" = "host",
           "OQ646790.1" = "phage", "OR222515.1" = "phage")

  p <- ggchord(seq_data_example, validate = "none") +
    geom_seq(
      seq_group = grp,
      seq_group_labels = c(host = "Host", phage = "Phage"),
      seq_group_colors = c(host = "firebrick", phage = "steelblue")
    )
  ggplot_build(p)
  layout <- get_chord_layout()

  expect_equal(layout$seq_groups, grp[layout$seqs])
  expect_equal(sort(layout$group_labels$label), c("Host", "Phage"))
  expected_colors <- unname(c("Host" = "firebrick", "Phage" = "steelblue")[
    layout$group_labels$label])
  expect_equal(unname(layout$group_labels$zcolour), expected_colors)
})

test_that("seq_group can be disabled and group labels can be hidden", {
  data(seq_data_example)
  seq_data_example$seq_group <- c("a", "a", "b", "b")

  p <- ggchord(seq_data_example, validate = "none") +
    geom_seq(seq_group_labels = FALSE)
  ggplot_build(p)
  layout <- get_chord_layout()
  expect_equal(nrow(layout$group_labels), 0)
  expect_equal(length(unique(layout$seq_groups)), 2)
})

test_that("seq_group errors clearly for invalid inputs", {
  data(seq_data_example)

  expect_error(
    ggplot_build(ggchord(seq_data_example, validate = "none") +
      geom_seq(seq_group = c("a", "b"))),
    "seq_group"
  )

  expect_error(
    ggplot_build(ggchord(seq_data_example, validate = "none") +
      geom_seq(seq_group_colors = c(not_a_group = "red"))),
    "seq_group_colors"
  )

  expect_error(
    ggplot_build(ggchord(seq_data_example, validate = "none") +
      geom_seq(seq_group_gap = -0.1)),
    "seq_group_gap"
  )
})

test_that("repeated non-contiguous groups get one label per contiguous block", {
  data(seq_data_example)
  seq_data_example$seq_group <- c("a", "b", "a", "b")

  p <- ggchord(seq_data_example, validate = "none") + geom_seq()
  ggplot_build(p)
  layout <- get_chord_layout()

  expect_equal(layout$group_labels$group_id, c("a", "b", "a", "b"))
  expect_equal(nrow(layout$group_labels), 4)
})
