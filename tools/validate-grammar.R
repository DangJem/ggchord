# Focused v0.16.1 public grammar checks. Run explicitly with Rscript.
# No testthat, package checks, rendered files or release device matrix.
pkgload::load_all(quiet = TRUE)
s <- data.frame(accver = c("A", "B"), length = c(1000, 1000))
g1 <- data.frame(accver = "A", start = 100, end = 300,
  strand = "+", anno = "one")
g2 <- data.frame(id = "A", lo = 500, hi = 700,
  direction = "-", category = "two")
features <- data.frame(accver = "A", start = 750, end = 900,
  strand = "+", type = "CDS", label = "feature")
results <- list()
record <- function(case, detail) {
  results[[length(results) + 1L]] <<- data.frame(case = case,
    result = "PASS", detail = detail)
}
for (mode in c("chord", "circular", "collinear")) {
  seqs <- if (mode == "circular") s[1, ] else s
  coord <- switch(mode, chord = coord_chord(), circular = coord_circular(),
    collinear = coord_collinear())
  scale <- scale_gene_fill_manual(values = c(one = "#113355", two = "#775533"))
  layers <- list(
    geom_gene(aes(gene_fill = anno), data = g1),
    geom_gene(aes(accver = id, start = lo, end = hi, strand = direction,
      anno = category, gene_fill = category), data = g2,
      position = position_strand(offset = .16)))
  p <- ggchord(seqs) + geom_seq() + scale + layers +
    geom_feature(data = features, shape = "block") +
    geom_gene_label(data = g1) + geom_feature_label(data = features) +
    geom_feature_label_repel(data = features) + coord
  b <- ggplot2::ggplot_build(p)
  stopifnot(nrow(b$data[[2]]) > 0, nrow(b$data[[3]]) > 0,
    all(b$data[[2]]$gene_fill == "#113355"),
    all(b$data[[3]]$gene_fill == "#775533"),
    all(vapply(b$data[4:7], nrow, integer(1)) > 0))
  record(paste0(mode, "_independent_gene_layers"),
    "Distinct data, role mappings, Position and manual colours retained")
  record(paste0(mode, "_feature_and_labels"),
    "Generic block feature, fixed gene/feature labels and feature repel build together")
  # A second build only changes layer/scale order, to detect shared layer state.
  reverse <- ggplot2::ggplot_build(ggchord(seqs) + geom_seq() +
    rev(layers) + scale + coord)
  stopifnot(all(reverse$data[[2]]$gene_fill == "#775533"),
    all(reverse$data[[3]]$gene_fill == "#113355"))
  record(paste0(mode, "_layer_and_scale_order"),
    "Reversing same-role layers and placing manual scale last retains each colour")
}
r <- data.frame(qaccver = "A", saccver = "B", qstart = c(100, 110),
  qend = c(200, 210), sstart = c(400, 410), send = c(500, 510),
  pident = c(90, 95))
for (mode in c("chord", "collinear")) {
  coord <- if (mode == "chord") coord_chord() else coord_collinear()
  p <- ggchord(s, r) + geom_seq() +
    stat_ribbon_density(aes(ribbon_alpha = after_stat(density),
      ribbon_fill = ggplot2::stage(start = pident,
        after_scale = rep("#224466", length(ribbon_fill)))),
      bins = 10, weight = "count") + coord
  b <- ggplot2::ggplot_build(p)
  stopifnot(nrow(b$data[[2]]) > 0,
    all(is.finite(b$data[[2]]$ribbon_alpha)),
    all(b$data[[2]]$ribbon_fill == "#224466"))
  record(paste0(mode, "_staged_density"),
    "after_stat density and staged after_scale evaluate on ribbon geometry")
}
sites <- data.frame(accver = "A", position = c(100, 120, 140, 450),
  enzyme = c("e1", "e2", "e3", "e4"))
geometry <- lapply(c("radial", "elbow", "straight"), function(style)
  get_chord_layout(ggchord(s[1, ]) + geom_seq() +
    geom_restriction_site(data = sites, leader = style) + coord_circular())[["restriction_sites"]])
leaders <- lapply(geometry, function(x) x[x$restriction_component == "leader", ])
vertices <- lapply(leaders, function(x) lengths(split(x$x, x$group)))
stopifnot(nrow(geometry[[1]]) > 0, all(vertices[[2]] == 3L),
  all(vertices[[3]] == 2L),
  identical(unique(leaders[[2]]$source_row), unique(leaders[[3]]$source_row)))
record("F39_leader_dispatch", "straight paths are direct; elbow paths bend once; source rows retained")
inside <- get_chord_layout(ggchord(s[1, ]) + geom_seq() +
  geom_restriction_site(data = sites, leader = "elbow", label_side = "inside") +
  coord_circular(direction = "counterclockwise", gap = 20))$restriction_sites
inside <- inside[inside$restriction_component == "leader", ]
paths <- split(inside, inside$group)
stopifnot(length(paths) > 0, all(vapply(paths, nrow, integer(1)) == 3L))
for (path in paths) stopifnot(sum(path[2, c("x", "y")]^2) <=
  sum(path[1, c("x", "y")]^2) + 1e-8)
record("F39_interior_elbow", "Interior elbow departs inward with reversed direction and a circular gap")
err <- tryCatch(theme_ggchord(axis.text.orientation = c(A = 0, B = 90)),
  error = identity)
stopifnot(inherits(err, "error"))
record("F40_public_orientation_domain", "Public theme rejects vectors; helper branch is outside public contract")
out <- do.call(rbind, results)
write.csv(out, "design-roadmap/audit-generated/grammar-validation-v0.16.1.csv",
  row.names = FALSE)
print(out, row.names = FALSE)
