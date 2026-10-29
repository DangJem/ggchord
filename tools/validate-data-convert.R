# Run explicitly from the repository root. No package test/check is performed.
# Set R_LIBS to an existing optional library. GGCHORD_REQUIRE_BIOC=true makes
# missing dependencies a failure rather than silently omitting their checks.
pkgload::load_all(quiet = TRUE)
local({
previous <- grDevices::dev.cur()
grDevices::pdf(NULL)
owned <- grDevices::dev.cur()
on.exit({
  if (owned %in% grDevices::dev.list()) grDevices::dev.off(owned)
  if (previous %in% grDevices::dev.list()) grDevices::dev.set(previous)
}, add = TRUE)
packages <- c("GenomicRanges", "InteractionSet", "Biostrings")
available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
if (identical(Sys.getenv("GGCHORD_REQUIRE_BIOC"), "true") && any(!available))
  stop("Required real-object validation dependencies unavailable: ",
    paste(packages[!available], collapse = ", "))

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
  exported <- export_ggchord_layout(direct,
    include = c("seq", "ribbon", "feature", "labels", "link"), original_data = TRUE)
  stopifnot(nrow(exported$feature) > 0L, nrow(exported$ribbon) > 0L,
    nrow(exported$link) > 0L, length(exported$conversion_reports) > 0L,
    exported$orientation_diagnostics$status[2L] == "no_evidence")
  grid::grid.draw(ggplot2::ggplotGrob(direct))
  unknown <- GenomicRanges::GRanges(seqnames = c("host", "virus"),
    ranges = IRanges::IRanges(start = c(100, 10), end = c(120, 30)))
  error <- tryCatch(as_ggchord_seq(unknown), error = identity)
  stopifnot(inherits(error, "error"),
    identical(as_ggchord_seq(unknown, lengths = c(host = 1000, virus = 300))$length,
      c(1000, 300)))
  cat("Real GRanges/GInteractions conversion checks passed.\n")
} else cat("Optional GRanges/GInteractions checks skipped: dependencies unavailable.\n")
if (available[["Biostrings"]]) {
  dna <- Biostrings::DNAStringSet(c("A.2" = paste(rep("ACGT", 25), collapse = ""),
    "B" = paste(rep("AT", 40), collapse = "")))
  S4Vectors::mcols(dna)$label <- c("host DNA", "viral DNA")
  S4Vectors::metadata(dna)$source <- "explicit synthetic adapter validation"
  converted <- as_ggchord_seq(dna)
  single <- as_ggchord_seq(Biostrings::DNAString("ACGT"), ids = "explicit")
  protein <- as_ggchord_seq(Biostrings::AAString("MKTW"), ids = "protein")
  characters <- as_ggchord_seq(Biostrings::BString("xY-1"), ids = "characters")
  stopifnot(identical(converted$accver, c("A.2", "B")),
    identical(converted$length, c(100, 80)), all(converted$sequence_unit == "bp"),
    identical(converted$label, c("host DNA", "viral DNA")),
    attr(converted, "ggchord_conversion")$source_metadata$source ==
      "explicit synthetic adapter validation",
    single$accver == "explicit", single$length == 4,
    protein$sequence_unit == "aa", protein$length == 4,
    characters$sequence_unit == "character", characters$length == 4)
  conflict <- tryCatch(as_ggchord_seq(dna, ids = c("renamed", "B")), error = identity)
  stopifnot(inherits(conflict, "error"))
  protein_plot <- tryCatch(ggchord(protein), error = identity)
  character_plot <- tryCatch(ggchord(characters), error = identity)
  stopifnot(inherits(protein_plot, "error"),
    grepl("project_ggchord_features", conditionMessage(protein_plot)),
    inherits(character_plot, "error"),
    grepl("known nucleotide lengths", conditionMessage(character_plot)))
  direct <- ggchord(dna) + geom_seq() + geom_seq_label() + coord_collinear()
  stopifnot(nrow(export_ggchord_layout(direct, include = "seq")$seq) > 0L)
  grid::grid.draw(ggplot2::ggplotGrob(direct))
  cat("Real XStringSet/XString names, widths, units, metadata and drawing checks passed.\n")
} else cat("Optional Biostrings checks skipped: dependency unavailable.\n")
for (name in packages[available]) cat(name, as.character(utils::packageVersion(name)), "\n")
cat("Data conversion development checks passed.\n")
})
