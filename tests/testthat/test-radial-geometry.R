count_radial_crossings <- function(paths) {
  if (nrow(paths) < 2L) return(0L)
  sum(vapply(seq_len(nrow(paths) - 1L), function(i) {
    j <- seq.int(i + 1L, nrow(paths))
    j <- j[paths$group[j] != paths$group[i]]
    sum(ggchord:::ggchord_segments_cross(paths$x0[i], paths$y0[i],
      paths$x1[i], paths$y1[i], paths$x0[j], paths$y0[j],
      paths$x1[j], paths$y1[j]))
  }, integer(1)))
}

test_that("nested measurements preserve the caller graphics device", {
  old_dir <- setwd(tempdir())
  on.exit(setwd(old_dir), add = TRUE)
  previous <- grDevices::dev.cur()
  original <- grDevices::dev.list()
  on.exit({
    for (device in setdiff(grDevices::dev.list(), original))
      grDevices::dev.off(device)
    if (previous %in% grDevices::dev.list()) grDevices::dev.set(previous)
  }, add = TRUE)
  grDevices::pdf(NULL, width = 6, height = 4)
  caller <- grDevices::dev.cur()
  p <- ggchord(smoke_sequences(), gene_data = smoke_genes(),
               validate = "none") +
    geom_seq() + geom_gene() + geom_gene_label_repel()
  before <- get_chord_layout(p)
  view_ggchord(p, viewer = "none")
  expect_identical(grDevices::dev.cur(), caller)
  expect_identical(get_chord_layout(p, build = FALSE)$gene_labels$label,
    before$gene_labels$label)
})


test_that("Unicode annotation measurement preserves content and caller device", {
  previous <- grDevices::dev.cur()
  annotation <- data.frame(text = "lacZ\u03b1", text_x = 0, text_y = 0, size = 2.5)
  boxes <- ggchord:::ggchord_text_boxes(annotation)
  expect_true(all(is.finite(boxes$w) & boxes$w > 0))
  expect_identical(grDevices::dev.cur(), previous)
})

test_that("radial labels retain inputs without overlapping or crossing", {
  seq <- data.frame(accver = "circle", length = 10000)
  genes <- data.frame(accver = "circle",
    start = seq(300, 9000, length.out = 6),
    end = seq(300, 9000, length.out = 6) + 100,
    strand = "+", anno = paste0("Site", seq_len(6)))
  p <- ggchord(seq, gene_data = genes, validate = "none") +
    geom_seq() + geom_gene() + geom_gene_label_repel()
  layout <- get_chord_layout(p)
  expect_equal(nrow(layout$gene_labels), nrow(genes))
  expect_false(ggchord:::ggchord_label_box_conflicts(layout$gene_labels,
    units_per_inch = layout$text_units_per_inch, box_padding = 0))
  expect_equal(count_radial_crossings(layout$gene_label_segments), 0L)
})

test_that("reversed orientation keeps leaders on their requested side", {
  seq <- data.frame(accver = c("A", "B"), length = 1000)
  genes <- data.frame(accver = rep(c("A", "B"), each = 2),
    start = c(200, 650, 200, 650), end = c(280, 730, 280, 730),
    strand = c("+", "-", "+", "-"), anno = letters[1:4])
  p <- ggchord(seq, gene_data = genes, validate = "none") +
    geom_seq(seq_orientation = c(1, -1)) + geom_gene() +
    geom_gene_label_repel(label_side = "outside")
  layout <- get_chord_layout(p)
  frame <- ggchord:::ggchord_label_curve_frame(layout$gene_labels,
                                                layout$seq_arcs)
  expect_true(all(frame$signed_distance > 0))
  expect_equal(count_radial_crossings(layout$gene_label_segments), 0L)
})

test_that("sparse restriction and primer leaders share circular root order", {
  s <- plasmid_example_pUC19
  primers <- find_primer_bindings(s, set = "universal")
  sites <- filter_restriction_sites(find_restriction_sites(s),
    enzymes = c("EcoRI", "BamHI", "HindIII"))
  p <- ggchord(s[c("accver", "length")]) + geom_seq() +
    geom_primer(data = primers) + geom_primer_label_repel(data = primers) +
    geom_restriction_site(data = sites) + coord_circular(rotation = 0)
  l <- get_chord_layout(p)
  paths <- list()
  for (id in names(l$layer_geometry)) {
    layer <- l$layer_geometry[[id]]
    g <- layer$gene_label_repel
    if (is.data.frame(g)) {
      g <- g[g$.component == "segment", ]
      if(nrow(g)) paths[[length(paths)+1L]] <- data.frame(x0=g$x,y0=g$y,
        x1=g$xend,y1=g$yend,group=paste(id,g$source_row))
    }
    g <- layer$restriction_site
    if (is.data.frame(g)) for (source in unique(g$source_row)) {
      z <- g[g$source_row == source & g$restriction_component == "leader", ]
      if(nrow(z)>1L) paths[[length(paths)+1L]] <- data.frame(
        x0=head(z$x,-1),y0=head(z$y,-1),x1=tail(z$x,-1),y1=tail(z$y,-1),
        group=paste(id,source))
    }
  }
  expect_equal(count_radial_crossings(do.call(rbind,paths)),0L)
})

test_that("explicit inside primer labels stay inside mixed circular maps", {
  s <- plasmid_example_pUC19
  primers <- find_primer_bindings(s, set = "universal")
  sites <- find_restriction_sites(s, enzymes = c("EcoRI", "BamHI", "HindIII"))
  p <- ggchord(s[c("accver", "length")]) + geom_seq() +
    geom_primer(data = primers) +
    geom_primer_label_repel(data = primers, side = "inside") +
    geom_restriction_site(data = sites) + coord_circular(rotation = 0)
  out <- export_ggchord_layout(p, include = "labels")
  text <- out$labels[out$labels$.component %in% "text", ]
  expect_setequal(text$source_row, seq_len(nrow(primers)))
  expect_true(all(sqrt(text$text_x^2 + text$text_y^2) <
    sqrt(text$anchor_x^2 + text$anchor_y^2)))
  registry <- out$annotation_registry
  expect_equal(sum(registry$kind == "primer_label" & registry$side == "inner"), nrow(primers))
})
