test_that("centre label lines stay in the panel with a reversed device y axis", {
  previous <- grDevices::dev.cur()
  grDevices::pdf(NULL, width = 7, height = 7)
  owned <- grDevices::dev.cur()
  on.exit({
    if (owned %in% grDevices::dev.list()) grDevices::dev.off(owned)
    if (previous %in% grDevices::dev.list()) grDevices::dev.set(previous)
  }, add = TRUE)
  plot <- ggchord(data.frame(accver = "circle", length = 1000),
    validate = "none") + geom_seq() + geom_seq_center_label() + coord_circular()
  # SVG and some screen devices reverse their native y axis. Reproduce that
  # coordinate convention on one private device without an export matrix.
  grid::pushViewport(grid::viewport(yscale = c(1, 0)))
  table <- tryCatch(ggplot2::ggplotGrob(plot), finally = grid::popViewport())
  collect <- function(grob) {
    if (inherits(grob, "text") &&
        any(grob$label %in% c("circle", "1,000 bp"))) return(list(grob))
    children <- if (inherits(grob, "gtable")) grob$grobs else grob$children
    unlist(lapply(children, collect), recursive = FALSE)
  }
  labels <- collect(table$grobs[[which(table$layout$name == "panel")]])
  grid::pushViewport(grid::viewport())
  y <- tryCatch(vapply(labels, function(grob) grid::convertY(
    grob$y, "npc", valueOnly = TRUE), numeric(1L)),
    finally = grid::popViewport())
  expect_length(y, 2L)
  expect_true(all(y > .25 & y < .75))
  expect_gt(y[1L], y[2L])
})


test_that("biological searches retain circular origin, strand and site identity", {
  primers <- find_primer_bindings(c(A = "TTGCAAAAACG"),
    primers = c(test = "ACGTTGCA"), circular = TRUE, min_annealed_bases = 8L)
  expect_true(primers$crosses_origin)
  expect_equal(primers$strand, "+")
  expect_error(find_primer_bindings(c(A = "ACGTACGTAAAATTTTGGGG"),
    primers = data.frame(primer_id = c("same", "same"),
      name = c("one", "two"), sequence = c("ACGTACGT", "TTTTGGGG"))),
    "each primer_id must identify one distinct sequence")
  db <- data.frame(common_feature_id = "gene", name = "gene", type = "CDS",
    directionality_label = "forward", sequence = "ATGAAA")
  hit <- find_common_features("TTTCAT", database = db, circular = FALSE)
  expect_equal(hit$strand, "-")
  sites <- find_restriction_sites(c(A = "TGAATTCG"), enzymes = "EcoRI")
  expect_equal(sites$enzyme, "EcoRI")
  expect_equal(sites$start, 2L)
})
