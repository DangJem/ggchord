test_that("v0.13 removes legacy sequence and ribbon entries", {
  old <- data.frame(seq_id = "A", length = 100)
  expect_error(ggchord(old), "removed in v0.13.0")
  expect_error(geom_seq(aes(seq_id = id)), "removed in v0.13.0")
  expect_false("geom_ribbon" %in% getNamespaceExports("ggchord"))
})

test_that("coord_circular owns a one-sequence circular contract", {
  seq <- data.frame(accver = "circle", length = 1000)
  expect_message(
    replaced <- ggchord(seq, validate = "none") + coord_circular(),
    NA
  )
  closed <- ggchord(seq, validate = "none") + geom_seq() + coord_circular()
  arc <- get_chord_layout(closed)$seq_arcs[[1L]]
  expect_s3_class(coord_circular(), "CoordCircular")
  expect_equal(
    unname(unlist(arc[1, c("x", "y")])),
    unname(unlist(arc[nrow(arc), c("x", "y")])), tolerance = 1e-8
  )

  open <- ggchord(seq, validate = "none") + geom_seq() +
    coord_circular(gap = 20, rotation = 90)
  open_arc <- get_chord_layout(open)$seq_arcs[[1L]]
  expect_gt(sum((open_arc[1, c("x", "y")] -
                   open_arc[nrow(open_arc), c("x", "y")])^2), .01)
  expect_equal(unname(unlist(open_arc[1, c("x", "y")])), c(0, 1), tolerance = .02)

  ccw <- ggchord(seq, validate = "none") + geom_seq() +
    coord_circular(gap = 20, direction = "counterclockwise")
  expect_identical(get_chord_layout(closed)$sequence_reference$orientation[[1]], -1)
  expect_identical(get_chord_layout(ccw)$sequence_reference$orientation[[1]], 1)
  expect_s3_class(ggplot2::ggplot_build(
    ggchord(seq, validate = "none") + geom_seq() + coord_circular(gap = 270)
  ), "ggplot_built")
  circular_build <- ggplot2::ggplot_build(closed)
  expect_null(circular_build$plot$layers[[1L]]$geom_params$arrow)
  expect_false(circular_build$plot$layers[[1L]]$show.legend)
  explicit_arrow <- ggplot2::ggplot_build(
    ggchord(seq, validate = "none") +
      geom_seq(arrow = grid::arrow(), show.legend = TRUE) + coord_circular()
  )
  expect_s3_class(explicit_arrow$plot$layers[[1L]]$geom_params$arrow, "arrow")
  expect_true(isTRUE(explicit_arrow$plot$layers[[1L]]$show.legend[["seq_colour"]]))

  two <- rbind(seq, transform(seq, accver = "other"))
  expect_error(
    get_chord_layout(ggchord(two, validate = "none") + geom_seq() + coord_circular()),
    "exactly one"
  )
  expect_error(
    get_chord_layout(ggchord(seq, validate = "none") +
      geom_seq(seq_gap = .1) + coord_circular()), "owns the opening"
  )
  expect_error(
    get_chord_layout(ggchord(seq, validate = "none") +
      geom_seq(seq_orientation = -1) + coord_circular()), "owns genomic direction"
  )
  ribbon <- data.frame(qaccver = "circle", saccver = "circle",
    qstart = 1, qend = 10, sstart = 20, send = 30)
  expect_error(
    get_chord_layout(ggchord(seq, ribbon, validate = "none") +
      geom_link_ribbon() + coord_circular()), "does not support"
  )
})

test_that("sequence backbone styles share the same reference", {
  seq <- data.frame(accver = "circle", length = 1000)
  build_style <- function(style) {
    p <- ggchord(seq, validate = "none") +
      geom_seq(seq_style = style) + coord_circular()
    ggplot2::ggplot_build(p)$data[[1L]]
  }
  single <- build_style("single")
  double <- build_style("double")
  band <- build_style("band")
  expect_true(all(single$.component == "path"))
  expect_equal(length(unique(double$group)), 2L)
  expect_true(all(double$.component == "path"))
  expect_true(all(band$.component == "band"))
  expect_error(geom_seq(seq_style = "double", seq_backbone_gap = -1), "backbone")
})

test_that("centre labels use sequence metadata", {
  seq <- data.frame(accver = "pBR322", label = "pBR322", length = 4361)
  p <- ggchord(seq, validate = "none") + geom_seq() +
    geom_seq_center_label() + coord_circular()
  center <- get_chord_layout(p)$seq_center_label
  expect_equal(center$label, "pBR322\n4,361 bp")
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
})

test_that("restriction search preserves biological pattern rows", {
  patterns <- data.frame(
    pattern_id = c("eco", "multi-a", "multi-b", "unknown", "type-iis", "four"),
    enzyme = c("EcoRI", "Multi", "Multi", "Mystery", "IIS", "FourCut"),
    motif = c("GAATTC", "RGATCY", "GGATCC", "AAAA", "GGTCTC", "CCGG"),
    ncuts = c(2L, 2L, 2L, 0L, 2L, 4L),
    blunt = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
    cut_offset_1 = c(1, 3, 1, 0, 7, -2),
    cut_offset_2 = c(5, 3, 5, 0, 11, -1),
    cut_offset_3 = c(0, 0, 0, 0, 0, 8),
    cut_offset_4 = c(0, 0, 0, 0, 0, 9),
    stringsAsFactors = FALSE
  )
  sites <- find_restriction_sites(
    c(g = "AAAAGAATTCTTTAGATCTGGATCCAAAAGGTCTCCCCGG"),
    patterns = patterns, circular = FALSE
  )
  expect_true(all(c(
    "match_id", "pattern_id", "pattern_source_row", "motif_length",
    "crosses_origin", "ncuts", "cut_offset_4", "cut_4_unwrapped",
    "cut_4", "display_position", "anchor_kind", "database_version"
  ) %in% names(sites)))
  expect_equal(length(unique(sites$pattern_id[sites$enzyme == "Multi"])), 2L)
  expect_true(any(sites$ncuts == 4L))
  expect_true(any(sites$ncuts == 0L & sites$anchor_kind == "recognition"))
  expect_true(any(sites$enzyme == "IIS" & sites$cut_2_unwrapped > sites$end))
  expect_true(any(sites$enzyme == "FourCut" & sites$cut_offset_1 < 0))
  eco <- sites[sites$enzyme == "EcoRI", , drop = FALSE]
  expect_equal(eco$cut_1_unwrapped, eco$start)
  expect_equal(eco$cut_2_unwrapped, eco$start + 4)

  directional <- data.frame(
    enzyme = "Directional", motif = "ACGTA", ncuts = 2L, blunt = FALSE,
    cut_offset_1 = 1L, cut_offset_2 = 4L
  )
  reverse_site <- find_restriction_sites(
    "TACGT", patterns = directional, circular = FALSE
  )
  expect_equal(reverse_site$strand, "-")
  expect_equal(reverse_site$cut_1_unwrapped, 1)
  expect_equal(reverse_site$cut_2_unwrapped, 4)

  wrap <- find_restriction_sites(
    c(circle = "AATTCG"), patterns = c(Wrap = "GAATTC"), circular = TRUE
  )
  expect_equal(nrow(wrap), 1L)
  expect_true(wrap$crosses_origin)
  expect_equal(wrap$start, 6L)
  expect_equal(nrow(find_restriction_sites(
    c(circle = "AATTCG"), patterns = c(Wrap = "GAATTC"), circular = FALSE
  )), 0L)

  duplicate <- rbind(
    transform(patterns[1, ], pattern_id = "dup-1"),
    transform(patterns[1, ], pattern_id = "dup-2")
  )
  dup_sites <- find_restriction_sites("GAATTC", patterns = duplicate)
  expect_equal(nrow(dup_sites), 2L)
  expect_equal(length(unique(dup_sites$pattern_id)), 2L)
})

test_that("REBASE parser uses stable pattern rows when source files exist", {
  path <- testthat::test_path("..", "..", "examples", "rebase")
  skip_if_not(all(file.exists(file.path(path, c(
    "VERSION", "embossa_e.txt", "embossa_r.txt", "embossa_s.txt"
  )))))
  parsed <- ggchord_parse_rebase(path)
  expect_equal(nrow(parsed), 4909L)
  expect_equal(length(unique(parsed$enzyme)), 4907L)
  expect_equal(sum(parsed$ncuts == 4L), 27L)
  expect_equal(sum(parsed$ncuts == 0L), 3264L)
  expect_true(any(parsed$cut_offset_1 < 0))
  expect_true(any(parsed$cut_offset_2 > parsed$motif_length))
  expect_true(all(grepl("^rebase609:e:[0-9]{6}$", parsed$pattern_id)))
})

test_that("restriction filters only subset rows", {
  sites <- data.frame(
    accver = "g", enzyme = c("A", "B", "B", "C"),
    motif_length = c(6L, 6L, 6L, 4L), position = c(10, 20, 30, 40),
    commercial = c(TRUE, TRUE, TRUE, FALSE), marker = seq_len(4)
  )
  expect_equal(filter_restriction_sites(sites, "unique")$enzyme, c("A", "C"))
  expect_equal(filter_restriction_sites(sites, "unique_dual")$marker, 1:4)
  expect_equal(filter_restriction_sites(sites, "six_plus")$marker, 1:3)
  expect_equal(filter_restriction_sites(sites, "unique_6plus")$marker, 1L)
  expect_equal(filter_restriction_sites(sites, "commercial")$marker, 1:3)
  out <- filter_restriction_sites(
    sites, set = "all", cuts = 2, window = c(15, 35), commercial_only = TRUE
  )
  expect_equal(out$marker, 2:3)
  expect_equal(out$position, sites$position[out$marker])
})

test_that("restriction layout is deterministic and never moves site anchors", {
  seq <- data.frame(accver = "g", length = 1000)
  sites <- data.frame(
    accver = "g", position = c(995, 5, 8, 700),
    enzyme = c("A", "B", "C", "D"), pattern_id = paste0("p", 1:4),
    pattern_source_row = 1:4
  )
  make_plot <- function(leader = "trunk") {
    ggchord(seq, validate = "none") + geom_seq() +
      geom_restriction_site(data = sites, leader = leader, min_label_gap = .02) +
      coord_circular(gap = 8)
  }
  first <- export_ggchord_layout(make_plot(), include = "restriction")$restriction
  second <- export_ggchord_layout(make_plot(), include = "restriction")$restriction
  expect_identical(first, second)
  ticks <- unique(first[first$restriction_component == "tick",
    c("source_row", "anchor_position")])
  expect_equal(ticks$anchor_position, sites$position[ticks$source_row])
  labels <- first[first$restriction_component == "label", , drop = FALSE]
  expect_equal(sort(unique(labels$source_row)), seq_len(nrow(sites)))
  expect_true(any(first$restriction_component == "trunk"))
  trunk <- first[first$restriction_component == "trunk", , drop = FALSE]
  expect_gt(nrow(trunk), 2L)
  expect_true(all(sqrt(trunk$x^2 + trunk$y^2) > 1))
  expect_s3_class(ggplot2::ggplot_build(make_plot("elbow")), "ggplot_built")
  expect_s3_class(ggplot2::ggplotGrob(make_plot("straight")), "gtable")
  expect_identical(
    export_ggchord_layout(make_plot(), include = "restriction")$metadata$coordinate,
    "circular"
  )
  styled <- ggplot2::ggplot_build(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_restriction_site(
        data = sites, colour = "black", linewidth = 1, label_size = 4
      ) + coord_circular() + theme_ggchord_plasmid()
  )$plot$layers[[2L]]$geom_params
  expect_equal(styled$segment_params$colour, "black")
  expect_equal(styled$segment_params$linewidth, 1)
  expect_equal(styled$text_params$size, 4)
})
