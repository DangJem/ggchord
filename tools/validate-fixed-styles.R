# Focused public behavior check; no package tests or graphics output.
pkgload::load_all(quiet = TRUE)
s <- data.frame(accver = c("A", "B"), length = c(1000, 1000))
r <- data.frame(qaccver = "A", saccver = "B", qstart = 100,
  qend = 200, sstart = 300, send = 400, pident = 90)
l <- data.frame(qaccver = "A", saccver = "B", qpos = 100, spos = 300,
  value = "mapped")
f <- data.frame(accver = "A", start = 600, end = 800,
  strand = "+", type = "CDS", shape_value = "block")
p <- ggchord(s, r) + geom_seq(aes(linewidth = I(.6)), colour = "#112233") +
  scale_link_colour_manual(values = c(mapped = "#112233")) +
  geom_link_line(aes(link_colour = value), data = l, colour = "#abcdef") +
  geom_link_line(aes(link_colour = value), data = l) +
  geom_link_ribbon(aes(ribbon_colour = "outline", ribbon_alpha = after_scale(.2)),
    colour = "#abcdef", alpha = .7) +
  geom_link_ribbon(aes(ribbon_colour = "outline", ribbon_alpha = pident)) +
  scale_ribbon_colour_manual(values = c(outline = "#112233")) +
  geom_feature(aes(feature_shape = shape_value), data = f,
    shape = "arrow") +
  geom_feature_plasmid(aes(feature_shape = shape_value), data = f) +
  geom_feature_plasmid(aes(feature_shape = shape_value), data = f,
    shape = "arrow") +
  geom_gene(aes(colour = I("#112233"), linewidth = I(.6)),
    data = transform(f, anno = "gene")) +
  geom_seq_region(aes(region_fill = I("#112233")), data = f,
    fill = "#abcdef", alpha = .8) + coord_collinear()
b <- ggplot2::ggplot_build(p)
stopifnot(all(b$data[[1]]$seq_colour == "#112233"),
  all(b$data[[1]]$linewidth == .6),
  all(b$data[[2]]$link_colour == "#abcdef"),
  all(b$data[[3]]$link_colour == "#112233"),
  all(b$data[[4]]$ribbon_colour == "#abcdef"),
  all(b$data[[4]]$ribbon_alpha == .7),
  all(b$data[[5]]$ribbon_colour == "#112233"),
  all(b$data[[6]]$feature_shape == "arrow"),
  all(b$data[[7]]$feature_shape == "block"),
  all(b$data[[8]]$feature_shape == "arrow"),
  all(b$data[[9]]$colour == "#112233"),
  all(b$data[[9]]$linewidth == .6),
  all(b$data[[10]]$region_fill == "#abcdef"),
  all(b$data[[10]]$alpha == .8))
inherited <- ggchord(s) + aes(link_colour = value) + geom_seq() +
  scale_link_colour_manual(values = c(mapped = "#112233")) +
  geom_link_line(data = l, inherit.aes = TRUE) +
  geom_link_line(data = l, inherit.aes = TRUE, colour = "#abcdef")
ib <- ggplot2::ggplot_build(inherited)
stopifnot(all(ib$data[[2]]$link_colour == "#112233"),
  all(ib$data[[3]]$link_colour == "#abcdef"))
cat("PASS: explicit fixed styles, omitted defaults, feature wrapper and inherited styles\n")
