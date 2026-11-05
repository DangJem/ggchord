# Scoped native-coordinate verification with observed protein inputs.
pkgload::load_all(quiet = TRUE)
proteins <- read_fasta_lengths(
  "examples/format-inputs/protein-pair.faa", sequence_unit = "aa"
)
hits <- read_blast("examples/format-inputs/protein-pair-blastp.tsv")
stopifnot(identical(proteins$length, c(296, 310)))
plain <- read_fasta_lengths("examples/format-inputs/protein-pair.faa")
stopifnot(!"sequence_unit" %in% names(plain))

p <- ggchord(proteins, hits) +
  geom_seq() +
  geom_link_ribbon() +
  coord_collinear() +
  geom_seq_scale(length = 50)
exported <- export_ggchord_layout(p, include = c("axis", "ribbon"))
stopifnot(
  exported$metadata$sequence_unit == "aa",
  any(exported$axis$label == "50 aa", na.rm = TRUE),
  all(is.na(exported$sequence_layout$units_per_bp)),
  all(is.finite(exported$sequence_layout$units_per_position))
)
aligned_interval <- data.frame(
  accver = hits$qaccver, start = hits$qstart, end = hits$qend,
  sequence_unit = "aa", directionality = "none"
)
invisible(ggplot2::ggplot_build(
  p + geom_feature(data = aligned_interval, shape = "block")
))
invisible(ggplot2::ggplot_build(
  ggchord(proteins, hits) + geom_seq() + geom_link_ribbon() + coord_chord()
))
circle <- ggchord(proteins[1, ]) +
  geom_seq() +
  geom_seq_center_label() +
  coord_circular()
text <- export_ggchord_layout(circle, include = "labels")$labels$label
stopifnot(any(grepl("296 aa", text)))

focused <- focus_ggchord_data(proteins, hits,
  loci = data.frame(accver = proteins$accver, start = 50, end = 200))
local <- ggchord(focused$seq_data, focused$ribbon_data) +
  geom_seq() +
  geom_link_ribbon() +
  coord_collinear()
point <- map_ggchord_position(local, proteins$accver[1], 100, source = TRUE)
back <- map_ggchord_position(local, proteins$accver[1], point$x,
  source = TRUE, inverse = TRUE)
stopifnot(back$position == 100, back$sequence_unit == "aa")
windows <- data.frame(accver = proteins$accver[1],
  start = c(20, 200), end = c(120, 280), locus_id = c("left", "right"))
pieces <- segment_ggchord_data(proteins[1, ], windows = windows)
invisible(ggplot2::ggplot_build(ggchord(pieces$seq_data) +
  geom_seq() + coord_collinear(segments = pieces$segments)))

expect_error_text <- function(expression, pattern) {
  result <- tryCatch(force(expression), error = identity)
  stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result)))
}
expect_error_text(ggchord(transform(proteins, sequence_unit = c("bp", "aa"))),
  "cannot share")
expect_error_text(ggchord(proteins, transform(hits, sequence_unit = "bp")),
  "disagrees")
expect_error_text(ggplot2::ggplot_build(p + geom_gene()), "require bp")
expect_error_text(find_restriction_sites(transform(proteins, sequence = "ACGT")),
  "DNA sequence input")
expect_error_text(find_primer_bindings(transform(proteins, sequence = "ACGT")),
  "DNA sequence input")
expect_error_text(find_common_features(transform(proteins, sequence = "ACGT")),
  "DNA sequence input")
if (requireNamespace("Biostrings", quietly = TRUE)) {
  aa <- Biostrings::readAAStringSet("examples/format-inputs/protein-pair.faa")
  stopifnot(get_chord_layout(ggchord(aa) + geom_seq())$sequence_unit == "aa")
}
cat("Native aa layouts, units, transformations and nucleotide-only boundaries passed.\n")
