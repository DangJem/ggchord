# One bounded public grob probe per layer/coordinate, using explicit synthetic geometry fixtures.
# This is documentation evidence, not a claim about all parameters or biological observations.
pkgload::load_all(quiet = TRUE)
sequences <- data.frame(accver = c("fixture_a", "fixture_b"), length = 1000)
features <- data.frame(accver = "fixture_a", start = 100, end = 300,
  strand = "+", anno = "fixture", feature_id = "fixture_feature")
links <- data.frame(qaccver = "fixture_a", saccver = "fixture_b",
  qstart = 100, qend = 300, sstart = 150, send = 350)
points <- data.frame(qaccver = "fixture_a", saccver = "fixture_b", qpos = 200, spos = 250)
primers <- data.frame(accver = "fixture_a", start = 120, end = 145,
  strand = "+", primer_id = "fixture_primer", name = "fixture_primer")
sites <- data.frame(accver = "fixture_a", position = 400, enzyme = "fixture_site")
bins <- data.frame(accver = "fixture_a", start = c(100, 300), end = c(199, 399), value = c(.2, .8))
samples <- data.frame(accver = "fixture_a", pos = c(150, 350), value = c(.2, .8))
constructors <- list(
  geom_seq = function() NULL,
  geom_seq_label = function() geom_seq_label(),
  geom_gene = function() geom_gene(),
  geom_feature = function() geom_feature(data = features),
  geom_feature_plasmid = function() geom_feature_plasmid(data = features),
  geom_gene_label = function() list(geom_gene(), geom_gene_label()),
  geom_gene_label_repel = function() list(geom_gene(), geom_gene_label_repel()),
  geom_feature_label = function() list(geom_feature(data = features), geom_feature_label(data = features)),
  geom_feature_label_repel = function() list(geom_feature(data = features), geom_feature_label_repel(data = features)),
  geom_link_ribbon = function() geom_link_ribbon(data = links),
  geom_link_line = function() geom_link_line(data = points),
  geom_seq_region = function() geom_seq_region(data = features),
  geom_seq_scale = function() geom_seq_scale(),
  geom_seq_break = function() geom_seq_break(),
  geom_seq_center_label = function() geom_seq_center_label(),
  geom_primer = function() geom_primer(data = primers),
  geom_primer_label_repel = function() list(geom_primer(data = primers), geom_primer_label_repel(data = primers)),
  geom_restriction_site = function() geom_restriction_site(data = sites),
  geom_track_bar = function() geom_track_bar(aes(track_value = value), data = bins, position = position_track("fixture")),
  geom_track_line = function() geom_track_line(aes(track_value = value), data = samples, position = position_track("fixture")),
  geom_track_tile = function() geom_track_tile(aes(track_fill = value), data = bins, position = position_track("fixture"))
)
results <- list()
for (mode in c("chord", "collinear", "circular")) {
  for (name in names(constructors)) {
    single <- mode == "circular" || name %in% c("geom_seq_center_label", "geom_primer", "geom_primer_label_repel", "geom_restriction_site")
    seqs <- if (single) sequences[1, , drop = FALSE] else sequences
    if (mode == "circular") {
      links$saccver <- "fixture_a"
      points$saccver <- "fixture_a"
    } else {
      links$saccver <- "fixture_b"
      points$saccver <- "fixture_b"
    }
    if (name == "geom_seq_break") {
      local <- focus_ggchord_data(seqs, loci = data.frame(accver = seqs$accver, start = 100, end = 900))
      seqs <- local$seq_data
    }
    coord <- switch(mode, chord = coord_chord(), collinear = coord_collinear(), circular = coord_circular())
    plot <- ggchord(seqs, gene_data = if (name == "geom_seq_break") NULL else features) + geom_seq() + coord
    result <- tryCatch({
      plot <- plot + constructors[[name]]()
      invisible(ggplot2::ggplotGrob(plot))
      list(status = "rendered", message = "")
    }, error = function(error) list(status = "rejected", message = conditionMessage(error)))
    results[[length(results) + 1L]] <- data.frame(api = name, coordinate = mode,
      sequences = nrow(seqs), status = result$status, message = result$message,
      fixture = "explicit synthetic geometry; no biological claims")
  }
}
write.csv(do.call(rbind, results), "tools/teaching-audit/layer-compatibility.csv", row.names = FALSE)
cat(length(results), "layer/coordinate grob outcomes recorded.\n")
