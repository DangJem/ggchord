# Shared finite RC inputs, not a new public API. Source after pkgload::load_all().
# No I/O or drawing at source time. Never alter bundled source records.
rc_fixtures <- function() {
  e <- new.env()
  for (name in c("seq_data_example", "ribbon_data_example", "gene_data_example",
    "plasmid_example_pUC19", "plasmid_example_pETDuet_1"))
    utils::data(list = name, package = "ggchord", envir = e)
  seqs <- e$seq_data_example[c("accver", "length")]
  genes <- transform(e$gene_data_example, strand = source_strand)
  genes <- genes[c("accver", "start", "end", "strand", "anno")]
  links <- e$ribbon_data_example[c("qaccver", "saccver", "qstart", "qend", "sstart", "send")]
  # Representative display points from real interval midpoints; not a new alignment.
  points <- data.frame(qaccver = links$qaccver, saccver = links$saccver,
    qpos = (links$qstart + links$qend) / 2, spos = (links$sstart + links$send) / 2)
  s <- e$plasmid_example_pUC19
  features <- find_common_features(s)
  features <- features[c("accver", "start", "end", "strand", "type", "anno")]
  primers <- find_primer_bindings(s, set = "universal")
  sites <- find_restriction_sites(s, enzymes = c("EcoRI", "BamHI", "HindIII"))
  # Bound the dense scene to true, nonredundant unique commercial cutters.
  # Selecting display rows never edits the original 49k+ recognition matches.
  dense_sites <- head(filter_restriction_sites(find_restriction_sites(s),
    set = "unique", parent_set = "commercial_nonredundant"), 40L)
  circular <- function(side = "outside", restriction_side = "outside", dense = FALSE) {
    ggchord(s[c("accver", "length")]) + geom_seq() +
      geom_feature(data = features, position = position_plasmid()) +
      geom_feature_label_repel(data = features, position = position_plasmid()) +
      geom_primer(data = primers) + geom_primer_label_repel(data = primers, side = side) +
      geom_restriction_site(data = if (dense) dense_sites else sites,
        label_side = restriction_side) +
      coord_circular(rotation = 90) + theme_ggchord_plasmid()
  }
  # Exact strand-corrected AmpR CDS equality supplies membership, not gene names.
  refs <- rbind(s, e$plasmid_example_pETDuet_1)
  cds <- do.call(rbind, lapply(seq_len(nrow(refs)), function(i) {
    f <- find_common_features(refs[i, ])
    f[f$type == "CDS" & f$anno == "AmpR", c("accver", "start", "end", "strand", "anno")]
  }))
  cds$feature_id <- with(cds, paste(accver, start, end, sep = ":"))
  dna <- vapply(seq_len(nrow(cds)), function(i) {
    f <- cds[i, ]; d <- substr(refs$sequence[match(f$accver, refs$accver)], f$start, f$end)
    if (f$strand == "-") d <- chartr("ACGT", "TGCA", paste(rev(strsplit(d, "")[[1L]]), collapse = ""))
    d
  }, character(1))
  stopifnot(nrow(cds) == 2L, length(unique(dna)) == 1L)
  homology <- link_ggchord_clusters(cds,
    data.frame(feature_id = cds$feature_id, cluster_id = "exact_CDS"),
    pairs = data.frame(qaccver = cds$accver[1L], saccver = cds$accver[2L]))$data
  windows <- transform(cds[c("accver", "start", "end")],
    start = pmax(1, start - 100), end = end + 100,
    locus_id = paste0("CDS", seq_len(nrow(cds))), row_id = accver,
    orientation = c(1, -1))
  view <- segment_ggchord_data(refs[c("accver", "length")], homology, cds,
    windows = windows)
  linear <- ggchord(view$seq_data, view$ribbon_data) + geom_seq() +
    geom_link_ribbon(fill = "#8971BB") + geom_feature(data = view$gene_data) +
    geom_feature_label_repel(data = view$gene_data) +
    coord_collinear(segments = view$segments) +
    scale_seq_position_continuous(coordinates = "source")
  starts <- seq.int(1L, s$length, by = 100L)
  bins <- data.frame(accver = s$accver, start = starts, end = pmin(starts + 99L, s$length))
  bins$value <- vapply(seq_len(nrow(bins)), function(i) {
    dna <- strsplit(substr(s$sequence, bins$start[i], bins$end[i]), "")[[1L]]
    mean(dna %in% c("G", "C"))
  }, numeric(1))
  track <- ggchord(s[c("accver", "length")]) + geom_seq() +
    geom_track_bar(aes(track_value = value, track_fill = value), data = bins,
      position = position_track("gc", width = unit(4, "mm"))) +
    scale_track_value_continuous("gc", limits = c(0, 1), unit = "GC fraction") +
    scale_track_fill_gradientn("gc", colours = c("#DCE9FF", "#214E95")) + coord_circular()
  list(seqs = seqs, genes = genes, links = links, points = points,
    source = s, features = features, primers = primers, sites = sites,
    circular = circular, cds = cds, homology = homology, windows = windows,
    view = view, linear = linear, bins = bins, track = track)
}

rc_leaders <- function(layout) {
  paths <- list()
  for (id in names(layout$layer_geometry)) {
    layer <- layout$layer_geometry[[id]]
    g <- layer$gene_label_repel
    if (is.data.frame(g)) {
      g <- g[g$.component %in% "segment", ]
      if (nrow(g)) paths[[length(paths) + 1L]] <- data.frame(x0 = g$x, y0 = g$y,
        x1 = g$xend, y1 = g$yend, group = paste(id, g$source_row))
    }
    g <- layer$restriction_site
    if (is.data.frame(g)) for (source in unique(g$source_row)) {
      z <- g[g$source_row == source & g$restriction_component %in% "leader", ]
      if (nrow(z) > 1L) paths[[length(paths) + 1L]] <- data.frame(
        x0 = head(z$x, -1), y0 = head(z$y, -1), x1 = tail(z$x, -1),
        y1 = tail(z$y, -1), group = paste(id, source))
    }
  }
  if (length(paths)) do.call(rbind, paths) else data.frame()
}

rc_crossings <- function(paths) {
  if (nrow(paths) < 2L) return(0L)
  sum(vapply(seq_len(nrow(paths) - 1L), function(i) {
    j <- seq.int(i + 1L, nrow(paths)); j <- j[paths$group[j] != paths$group[i]]
    sum(ggchord:::ggchord_segments_cross(paths$x0[i], paths$y0[i], paths$x1[i],
      paths$y1[i], paths$x0[j], paths$y0[j], paths$x1[j], paths$y1[j]))
  }, integer(1)))
}
