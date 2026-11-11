# Targeted development checks, not package-wide release acceptance.
# Run from repository root: Rscript tools/long-term-audit/check-dotplot.R
pkgload::load_all(quiet = TRUE)
library(ggplot2)

check_endpoints <- function(seqs, hits, name) {
  original <- hits
  p <- ggchord(seqs, hits, validate = "none") +
    geom_alignment_segment(aes(link_colour = pident)) + coord_dotplot()
  exported <- export_ggchord_layout(p, original_data = TRUE)
  d <- exported$alignment
  stopifnot(nrow(d) == nrow(hits), identical(hits, original))
  for (axis in c("x", "y")) {
    prefix <- if (axis == "x") "q" else "s"
    for (end in c("start", "end")) {
      display <- if (end == "start") axis else paste0(axis, "end")
      back <- map_ggchord_position(p, d[[paste0(prefix, "accver")]],
        d[[display]], axis = axis, inverse = TRUE)
      stopifnot(isTRUE(all.equal(back$position, d[[paste0(prefix, end)]])))
    }
  }
  built <- ggplot_build(p)
  stopifnot(sum(vapply(built$data, nrow, integer(1))) >= nrow(hits))
  invisible(ggplot_gtable(built))
  message(name, ": ", nrow(hits), " HSP endpoints and rendered build passed")
  p
}

# Three real source-backed cases; PAF declared lengths, nucleotide FASTA,
# and protein FASTA respectively. No lengths or identity values are invented.
paf <- read_paf("examples/format-inputs/emales-pair.paf")
seqs <- unique(rbind(data.frame(accver = paf$qaccver, length = paf$qlen),
  data.frame(accver = paf$saccver, length = paf$slen)))
real <- check_endpoints(seqs, paf, "EMALE PAF")
nucl <- read_blast("examples/blastn-standard/MT108731.1__MT118296.1.tsv")
seqs <- read_fasta_lengths(files = c("examples/fasta/MT108731.1.fna", "examples/fasta/MT118296.1.fna"))
invisible(check_endpoints(seqs, nucl, "Nucleotide BLAST"))
protein <- read_blast("examples/format-inputs/protein-pair-blastp.tsv")
seqs <- read_fasta_lengths("examples/format-inputs/protein-pair.faa")
seqs$sequence_unit <- "aa"; protein$sequence_unit <- "aa"
invisible(check_endpoints(seqs, protein, "Protein BLAST"))

# Explicitly synthetic boundary case: independent axis order/orientation,
# two independent layers with different metadata, duplicate HSP and self hit.
seqs <- data.frame(accver = c("A", "B", "C"), length = c(1000, 800, 300))
hits <- data.frame(qaccver = c("A", "B", "A"), saccver = c("B", "A", "A"),
  qstart = c(1, 10, 40), qend = c(500, 700, 80), sstart = c(800, 20, 40),
  send = c(100, 900, 80), pident = c(95, 80, 100), source_hsp_id = c("h1", "h2", "self"))
extra <- hits[1, ]; extra$note <- "independent second layer"
p <- ggchord(seqs, hits, validate = "none") +
  geom_alignment_segment(aes(link_colour = pident), allow_repeated = TRUE) +
  geom_alignment_segment(data = extra, colour = "red", allow_repeated = TRUE) +
  coord_dotplot(x_seqs = c("B", "A", "C"), y_seqs = c("A", "B"),
    x_orientation = c(B = -1, A = 1, C = 1), y_orientation = -1, gap = 100)
x <- export_ggchord_layout(p, original_data = TRUE)
stopifnot(nrow(x$alignment) == 4L, length(unique(x$alignment$layer_id)) == 2L,
  identical(x$original_data[[1]]$alignment, hits))
for (axis in c("x", "y")) {
  d <- x$alignment
  side <- if (axis == "x") "q" else "s"
  back <- map_ggchord_position(p, d[[paste0(side, "accver")]], d[[axis]], axis = axis, inverse = TRUE)
  stopifnot(isTRUE(all.equal(back$position, d[[paste0(side, "start")]])))
}
stopifnot(is.na(map_ggchord_position(p, "B", 850, inverse = TRUE, axis = "x")$position))
# A user-supplied style scale wins, and native theme overrides remain effective.
custom <- p + scale_link_colour_gradientn(colours = c("black", "white"), limits = c(70, 100)) +
  theme(axis.text.x = element_blank())
b <- ggplot_build(custom)
stopifnot(identical(b$plot$scales$get_scales("link_colour")$limits, c(70, 100)),
  inherits(calc_element("axis.text.x", b$plot$theme), "element_blank"))
invisible(ggplot_gtable(b))
# Selection omissions, unknown direction and contradictory metadata.
selected <- export_ggchord_layout(ggchord(seqs, hits, validate = "none") +
  geom_alignment_segment() + coord_dotplot(x_seqs = "B", y_seqs = "A"))
stopifnot(sum(selected$connection_report$retained) == 1L)
unknown <- hits[1, ]; unknown$direction <- "unknown"
u <- suppressWarnings(export_ggchord_layout(ggchord(seqs, unknown, validate = "none") +
  geom_alignment_segment() + coord_dotplot()))
stopifnot(nrow(u$alignment) == 0L, u$connection_report$reason == "unknown_direction")
unknown$direction <- "same"
err <- tryCatch({export_ggchord_layout(ggchord(seqs, unknown, validate = "none") +
  geom_alignment_segment() + coord_dotplot()); FALSE}, error = function(e) grepl("direction conflicts", conditionMessage(e)))
stopifnot(err)
# Focused source origins remain reversible on independently reversed axes.
focused <- focus_ggchord_data(seqs, hits, loci = data.frame(
  accver = c("A", "B"), start = c(20, 10), end = c(950, 790)))
fp <- ggchord(focused$seq_data, focused$ribbon_data, validate = "none") +
  geom_alignment_segment() + coord_dotplot(x_orientation = -1)
fpos <- map_ggchord_position(fp, "A", 100, axis = "x", source = TRUE)
stopifnot(map_ggchord_position(fp, "A", fpos$x, axis = "x", source = TRUE, inverse = TRUE)$position == 100)
# Default requires explicit authorization when a declared HSP is repeated.
repeated <- tryCatch({get_chord_layout(ggchord(seqs, hits, validate = "none") +
  geom_alignment_segment() + geom_alignment_segment(data = extra) + coord_dotplot()); FALSE},
  error = function(e) grepl("allow_repeated", conditionMessage(e)))
stopifnot(repeated)
# Empty alignment layers still have complete selected sequence axes.
empty <- ggchord(seqs) + geom_alignment_segment() + coord_dotplot(x_seqs = "A", y_seqs = "B")
stopifnot(nrow(export_ggchord_layout(empty)$alignment) == 0L)
invisible(ggplot_gtable(ggplot_build(empty)))
# Existing view still builds after shared entry-point changes.
invisible(ggplot_build(ggchord(seqs, hits, validate = "none") + geom_seq() + geom_link_ribbon()))
message("Synthetic axis/provenance/style/omission checks and existing Chord build passed")
if (identical(Sys.getenv("GGCHORD_DOTPLOT_PREVIEW"), "1"))
  ggsave(file.path(tempdir(), "ggchord-dotplot-real.png"), real, width = 7, height = 6)
