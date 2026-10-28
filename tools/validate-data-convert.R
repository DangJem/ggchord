# Run explicitly from the repository root. No package test/check is performed.
pkgload::load_all(quiet = TRUE)

seqs <- as_ggchord_seq(data.frame(chrom = c("host", "virus"), size = c(1000, 300)),
  mapping = c(accver = "chrom", length = "size"))
stopifnot(identical(seqs$accver, c("host", "virus")), all(seqs$length > 0))
features <- as_ggchord_feature(data.frame(chrom = "host", from = 100, to = 200,
  direction = "*", label = "integration", .source_row = 99L),
  mapping = c(accver = "chrom", start = "from", end = "to", strand = "direction"))
stopifnot(features$strand == ".", features$directionality == "none",
  features$source_strand == "*", features$metadata_.source_row == 99L)
links <- as_ggchord_link(data.frame(host = "host", virus = "virus",
  host_start = 100, host_end = 120, virus_start = 30, virus_end = 10,
  score = 12, qaccver = "source metadata"),
  mapping = c(qaccver = "host", saccver = "virus", qstart = "host_start",
    qend = "host_end", sstart = "virus_start", send = "virus_end"))
stopifnot(links$qaccver == "host", links$metadata_qaccver == "source metadata",
  links$sstart > links$send, links$score == 12,
  !any(c("pident", "length") %in% names(links)),
  identical(attr(links, "ggchord_conversion")$coordinate_system, "1-based-inclusive"))
plot <- ggchord(seqs, links, validate = "none") + geom_seq() + geom_link_ribbon() +
  geom_feature(data = features) + coord_collinear()
out <- export_ggchord_layout(plot, include = c("seq", "ribbon", "feature"))
stopifnot(nrow(out$feature) > 0L, nrow(out$ribbon) > 0L)
bad_length <- tryCatch(as_ggchord_seq(data.frame(accver = "host")), error = identity)
stopifnot(inherits(bad_length, "error"), grepl("full sequence lengths", conditionMessage(bad_length)))
bad_origin <- tryCatch(as_ggchord_feature(data.frame(accver = "host", start = 0, end = 10)), error = identity)
stopifnot(inherits(bad_origin, "error"))

if (requireNamespace("GenomicRanges", quietly = TRUE) &&
    requireNamespace("InteractionSet", quietly = TRUE)) {
  gr <- GenomicRanges::GRanges(seqnames = c("host", "virus"),
    ranges = IRanges::IRanges(start = c(100, 10), end = c(120, 30)),
    strand = c("-", "*"), seqlengths = c(host = 1000L, virus = 300L))
  names(gr) <- c("host interval", "virus interval")
  S4Vectors::mcols(gr)$label <- c("host gene", "viral site")
  S4Vectors::mcols(gr)$accver <- c("metadata host", "metadata virus")
  S4Vectors::metadata(gr)$reference <- "synthetic validation intervals"
  gf <- as_ggchord_feature(gr)
  gs <- as_ggchord_seq(gr)
  info_seq <- as_ggchord_seq(GenomicRanges::seqinfo(gr))
  stopifnot(identical(gf$accver, c("host", "virus")), gf$strand[2L] == ".",
    gf$directionality[2L] == "none", gf$metadata_accver[1L] == "metadata host",
    identical(gs$length, c(1000, 300)), identical(info_seq$accver, gs$accver))
  gi <- InteractionSet::GInteractions(gr[2L], gr[1L])
  S4Vectors::mcols(gi)$score <- 0.8
  S4Vectors::mcols(gi)$qaccver <- "interaction metadata"
  gl <- as_ggchord_link(gi)
  stopifnot(gl$qaccver == "virus", gl$saccver == "host", gl$sstart == 120,
    gl$send == 100, !gl$orientation_known, gl$score == .8,
    gl$metadata_qaccver == "interaction metadata",
    !any(c("pident", "length") %in% names(gl)))
  # Optional S4 inputs also work directly through the constructor and layers.
  direct <- ggchord(GenomicRanges::seqinfo(gr), gi, validate = "none") +
    geom_seq() + geom_link_ribbon() + geom_feature(data = gr) +
    geom_feature_label_repel(data = gr) + geom_link_line(data = gi) +
    coord_collinear(orientation = "auto")
  exported <- export_ggchord_layout(direct, original_data = TRUE)
  stopifnot(nrow(exported$feature) > 0L, nrow(exported$ribbon) > 0L,
    nrow(exported$link) > 0L, length(exported$conversion_reports) > 0L,
    exported$orientation_diagnostics$status[2L] == "no_evidence")
  unknown <- GenomicRanges::GRanges(seqnames = c("host", "virus"),
    ranges = IRanges::IRanges(start = c(100, 10), end = c(120, 30)))
  error <- tryCatch(as_ggchord_seq(unknown), error = identity)
  stopifnot(inherits(error, "error"),
    identical(as_ggchord_seq(unknown, lengths = c(host = 1000, virus = 300))$length,
      c(1000, 300)))
  cat("Real GRanges/GInteractions conversion checks passed.\n")
} else cat("Optional GRanges/GInteractions checks skipped: dependencies unavailable.\n")
cat("Data conversion development checks passed.\n")
