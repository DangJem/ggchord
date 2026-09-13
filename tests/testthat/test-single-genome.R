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

  plasmid_axis <- get_chord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      scale_seq_position_continuous() + coord_circular() +
      theme_ggchord_plasmid()
  )
  major <- plasmid_axis$axis_ticks[plasmid_axis$axis_ticks$is_major, ]
  expect_lt(
    mean(sqrt(major$label_x^2 + major$label_y^2)),
    mean(sqrt(plasmid_axis$seq_arcs[[1L]]$x^2 +
      plasmid_axis$seq_arcs[[1L]]$y^2))
  )

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
  auto <- build_style("auto")
  double <- build_style("double")
  band <- build_style("band")
  expect_true(all(single$.component == "path"))
  expect_equal(length(unique(auto$group)), 2L)
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
  expect_equal(center$label, c("pBR322", "4,361 bp"))
  expect_equal(center$.component, c("name", "length"))
  expect_true(all(center$combined_label == "pBR322\n4,361 bp"))
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")

  styled <- ggchord(seq, validate = "none") + geom_seq() +
    geom_seq_center_label(
      name_style = list(colour = "red", fontface = "bold"),
      length_style = list(size = 2.5)
    ) + coord_circular()
  expect_s3_class(ggplot2::ggplotGrob(styled), "gtable")
})

test_that("external feature and restriction labels share perimeter space", {
  seq <- data.frame(accver = "circle", length = 1000)
  features <- data.frame(
    accver = "circle", start = c(95, 120, 145), end = c(101, 127, 153),
    strand = "+", anno = paste("short feature", 1:3),
    feature_color = c("#31849B", "#FF0000", "#FFFFFF")
  )
  sites <- data.frame(
    accver = "circle", position = c(90, 135, 160, 650),
    enzyme = c("OnceA", "OnceB", "Repeat", "Repeat")
  )
  tracks <- position_feature_stack(base_position = position_plasmid())
  plot <- ggchord(seq, validate = "none") + geom_seq() +
    geom_feature_plasmid(data = features, position = tracks) +
    geom_feature_label_repel(data = features, position = tracks) +
    geom_restriction_site(data = sites) + coord_circular(rotation = 90)
  exported <- export_ggchord_layout(
    plot, include = c("labels", "restriction")
  )
  feature <- exported$labels[
    exported$labels$.component == "text" &
      exported$labels$feature_label_mode == "external", , drop = FALSE
  ]
  restriction <- exported$restriction[
    exported$restriction$.component == "label", , drop = FALSE
  ]
  expect_true(nrow(feature) > 0L)
  expect_true(all(feature$external_annotation_type == "feature"))
  expect_true(all(restriction$external_annotation_type == "restriction"))
  expect_true(all(restriction$enzyme_fontface[
    restriction$enzyme_label %in% c("OnceA", "OnceB")
  ] == "bold"))
  expect_true(all(restriction$enzyme_fontface[
    restriction$enzyme_label == "Repeat"
  ] == "plain"))

  external <- rbind(
    data.frame(text_x = feature$x, text_y = feature$y,
      text = feature$label, size = feature$size, hjust = feature$hjust,
      vjust = feature$vjust, text_angle = 0),
    data.frame(text_x = restriction$x, text_y = restriction$y,
      text = restriction$label, size = restriction$size,
      hjust = restriction$hjust, vjust = restriction$vjust,
      text_angle = 0)
  )
  boxes <- ggchord_text_boxes(
    external, units_per_inch = get_chord_layout(plot)$text_units_per_inch,
    box_padding = .01
  )
  if (nrow(boxes) > 1L) {
    pairs <- utils::combn(seq_len(nrow(boxes)), 2L)
    expect_false(any(apply(pairs, 2L, function(pair) {
      ggchord_oriented_box_overlaps(
        boxes[pair[1L], , drop = FALSE], boxes[pair[2L], , drop = FALSE]
      )
    })))
  }
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
    "cut_4", "display_position", "enzyme_site_count", "anchor_kind",
    "database_version"
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

  overlapping <- find_restriction_sites(
    "GGGG", patterns = c(Overlap = "GGG"), circular = FALSE
  )
  expect_equal(overlapping$start, 1:2)
  expect_equal(unique(overlapping$enzyme_site_count), 2L)
  one_visible <- filter_restriction_sites(overlapping, window = c(1, 1))
  expect_equal(nrow(one_visible), 1L)
  expect_equal(one_visible$enzyme_site_count, 2L)
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
  expect_true(all(c("source_motif", "preferred_enzyme",
    "is_preferred_enzyme") %in% names(parsed)))
  expect_identical(parsed$preferred_enzyme[parsed$enzyme == "PspFI"], "BseYI")
  expect_true(any(parsed$source_motif != parsed$motif))
  expect_equal(parsed$cut_offset_1[parsed$enzyme == "PspFI"], 5L)
  expect_equal(parsed$cut_offset_2[parsed$enzyme == "PspFI"], 1L)
  expect_true(all(grepl("^rebase609:e:[0-9]{6}$", parsed$pattern_id)))
})

test_that("pBluescript IUPAC and reverse-cleavage sites match references", {
  data(plasmid_example_pBluescript_II_SK_plus)
  sites <- find_restriction_sites(plasmid_example_pBluescript_II_SK_plus,
    enzymes = c("EcoO109I", "PspFI", "BseYI"))
  observed <- sites$position[match(c("EcoO109I", "PspFI", "BseYI"),
    sites$enzyme)]
  expect_equal(observed, c(660, 1461, 1457))
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

  equivalent <- data.frame(
    accver = "g", enzyme = c("Alias", "Preferred", "Other"),
    preferred_enzyme = c("Preferred", "Preferred", "Other"),
    is_preferred_enzyme = c(FALSE, TRUE, TRUE),
    pattern_source_row = 1:3, motif_length = 6L, start = c(10L, 10L, 20L),
    position = c(10, 11, 20), commercial = TRUE
  )
  reduced <- filter_restriction_sites(equivalent,
    parent_set = "commercial_nonredundant")
  expect_equal(reduced$enzyme, c("Preferred", "Other"))
})

test_that("restriction layout is deterministic and never moves site anchors", {
  seq <- data.frame(accver = "g", length = 1000)
  sites <- data.frame(
    accver = "g", position = c(995, 5, 8, 700),
    enzyme = c("A", "B", "C", "D"), pattern_id = paste0("p", 1:4),
    pattern_source_row = 1:4
  )
  make_plot <- function(leader = "radial") {
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
  expect_false(any(first$restriction_component == "trunk"))
  leaders <- first[first$restriction_component == "leader", , drop = FALSE]
  expect_true(nrow(leaders) > 0L)
  expect_true(all(lengths(leaders$source_rows) == 1L))
  expect_true(all(c(
    "source_rows", "cluster_id", "junction_id", "label_direction",
    "label_order", "label_connection_side", "plotmath_label"
  ) %in% names(first)))
  expect_true(all(labels$label_order[
    labels$label_connection_side == "left"] == "enzyme_position"))
  expect_true(all(labels$label_order[
    labels$label_connection_side == "right"] == "position_enzyme"))
  # Validate the final rendered approach, not just hjust: a right-travelling
  # segment enters the left edge and a left-travelling segment the right edge.
  for (junction in labels$junction_id) {
    label_row <- labels[labels$junction_id == junction, , drop = FALSE]
    leader_rows <- leaders[leaders$junction_id == junction, , drop = FALSE]
    # Each segment is exported as its own ordered path; the final path has the
    # largest group id. A bend may be closer to a long label's centre than its
    # true edge endpoint, so nearest-point selection is not sufficient here.
    final_group <- max(leader_rows$group)
    final_path <- leader_rows[leader_rows$group == final_group, ]
    distance <- (final_path$x - label_row$x)^2 +
      (final_path$y - label_row$y)^2
    endpoint <- final_path[which.min(distance), , drop = FALSE]
    start <- final_path[which.max(
      (final_path$x - label_row$x)^2 + (final_path$y - label_row$y)^2
    ), , drop = FALSE]
    if (label_row$label_connection_side == "left") {
      expect_lt(endpoint$x, label_row$x)
      expect_gt(endpoint$x, start$x)
    } else {
      expect_gt(endpoint$x, label_row$x)
      expect_lt(endpoint$x, start$x)
    }
  }
  expect_true(all(labels$label_connection_side %in% c("left", "right")))
  expect_s3_class(ggplot2::ggplot_build(make_plot("elbow")), "ggplot_built")
  expect_s3_class(ggplot2::ggplotGrob(make_plot("straight")), "gtable")
  ordering <- function(style) {
    x <- export_ggchord_layout(make_plot(style), include = "restriction")$restriction
    x <- x[x$restriction_component == "label",
      c("source_row", "slot_position", "label_direction", "label_order")]
    x <- x[order(x$source_row), , drop = FALSE]
    rownames(x) <- NULL
    x
  }
  expect_identical(ordering("radial"), ordering("elbow"))
  expect_identical(ordering("radial"), ordering("straight"))
  expect_identical(
    export_ggchord_layout(make_plot(), include = "restriction")$metadata$coordinate,
    "circular"
  )
  origin_labels <- labels[labels$anchor_position %in% c(995, 5, 8), ]
  origin_labels <- origin_labels[order(-origin_labels$y), ]
  expect_equal(origin_labels$anchor_position, c(995, 5, 8))
  origin_radius <- sqrt(origin_labels$x^2 + origin_labels$y^2)
  expect_lt(diff(range(origin_radius)), .04)

  sparse <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_restriction_site(data = sites[4, , drop = FALSE]) +
      coord_circular(), include = "restriction"
  )$restriction
  sparse <- sparse[sparse$restriction_component == "leader", , drop = FALSE]
  expect_lte(length(unique(sparse$group)), 2L)

  dense_sites <- data.frame(
    accver = "g", position = seq(1, 22, by = 3),
    enzyme = paste0("LongEnzyme", seq_len(8))
  )
  dense <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_restriction_site(data = dense_sites) + coord_circular(),
    include = "restriction"
  )$restriction
  dense_leaders <- dense[dense$restriction_component == "leader", ]
  segments_per_site <- vapply(split(dense_leaders, dense_leaders$source_row),
    function(x) length(unique(x$group)), integer(1))
  expect_true(all(segments_per_site <= 2L))
  expect_true(any(segments_per_site == 2L))

  lateral_sites <- data.frame(
    accver = "g", position = seq(220, 248, by = 4),
    enzyme = paste0("E", seq_len(8))
  )
  lateral <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_restriction_site(data = lateral_sites) +
      coord_circular(rotation = 90),
    include = "restriction"
  )$restriction
  lateral_labels <- lateral[
    lateral$restriction_component == "label", , drop = FALSE
  ]
  expect_true(all(lateral_labels$label_layout == "contour_fan"))
  boundary_radius <- sqrt(
    lateral_labels$label_boundary^2 + lateral_labels$y^2
  )
  expect_equal(boundary_radius, rep(boundary_radius[1L],
    nrow(lateral_labels)), tolerance = 1e-8)
  expect_gt(length(unique(round(lateral_labels$label_boundary, 4))), 2L)
  ordered_column <- lateral_labels[order(lateral_labels$anchor_position), ]
  expect_true(all(diff(ordered_column$y) < 0))
  expect_equal(abs(diff(ordered_column$y)),
    rep(abs(diff(ordered_column$y))[1L], nrow(ordered_column) - 1L),
    tolerance = 1e-8)
  expect_gt(min(abs(diff(ordered_column$y))), .045)
  column_leaders <- lateral[
    lateral$restriction_component == "leader", , drop = FALSE
  ]
  terminal_x <- bend_x <- numeric(nrow(lateral_labels))
  for (i in seq_len(nrow(lateral_labels))) {
    paths <- column_leaders[
      column_leaders$junction_id == lateral_labels$junction_id[i], ]
    final <- paths[paths$group == max(paths$group), ]
    terminal_x[i] <- final$x[which.min(
      (final$x - lateral_labels$x[i])^2 +
        (final$y - lateral_labels$y[i])^2
    )]
    point_key <- paste(round(paths$x, 10), round(paths$y, 10), sep = "\r")
    repeated <- paths$x[duplicated(point_key)]
    bend_x[i] <- if (length(repeated)) repeated[1L] else NA_real_
  }
  expect_equal(terminal_x, lateral_labels$label_boundary - .003,
    tolerance = 1e-8)
  expect_gte(sum(is.finite(bend_x)), 2L)
  expect_gt(diff(range(bend_x, na.rm = TRUE)), .005)

  small_cluster <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_restriction_site(data = lateral_sites[1:4, ]) +
      coord_circular(rotation = 90),
    include = "restriction"
  )$restriction
  small_cluster <- small_cluster[
    small_cluster$restriction_component == "label", , drop = FALSE
  ]
  expect_true(all(small_cluster$label_layout == "contour"))
  attachment_radius <- sqrt(
    small_cluster$label_attachment_x^2 +
      small_cluster$label_attachment_y^2
  )
  expect_equal(attachment_radius, small_cluster$label_contour_radius,
    tolerance = 1e-8)
  expect_lt(diff(range(attachment_radius)), 1e-8)

  same_position <- data.frame(
    accver = "g", position = c(250, 250),
    enzyme = c("BsaAI", "DraIII"),
    pattern_id = c("p1", "p2"), pattern_source_row = 1:2
  )
  combined <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_restriction_site(data = same_position) + coord_circular(),
    include = "restriction"
  )$restriction
  combined_label <- combined[
    combined$restriction_component == "label", , drop = FALSE
  ]
  expect_equal(nrow(combined_label), 1L)
  expect_match(combined_label$label, "BsaAI - DraIII")
  expect_equal(combined_label$enzyme_label, "BsaAI - DraIII")
  expect_equal(combined_label$coordinate_label, "(250)")
  expect_equal(combined_label$source_rows[[1L]], 1:2)
  expect_equal(length(unique(combined$group[
    combined$restriction_component == "tick"
  ])), 1L)

  neighbouring <- same_position
  neighbouring$position <- c(250, 251)
  separate <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_restriction_site(data = neighbouring) + coord_circular(),
    include = "restriction"
  )$restriction
  expect_equal(sum(separate$restriction_component == "label"), 2L)

  perimeter_sites <- data.frame(
    accver = "g", position = seq(25, 975, length.out = 16),
    enzyme = paste0("E", seq_len(16))
  )
  perimeter <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_restriction_site(data = perimeter_sites) + coord_circular(),
    include = "restriction"
  )$restriction
  perimeter <- perimeter[perimeter$restriction_component == "label", ]
  for (direction in c("left", "right")) {
    side <- perimeter[perimeter$label_direction == direction, ]
    expect_gt(length(unique(round(side$x, 4))), 2L)
  }
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
