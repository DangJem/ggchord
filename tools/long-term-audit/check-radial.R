# Targeted Radial development checks. No package-wide release checks.
# Run from repo root: Rscript tools/long-term-audit/check-radial.R
pkgload::load_all(quiet = TRUE)
library(ggplot2)
fail <- function(expr, pattern) {
  message <- tryCatch({force(expr); "NO ERROR"}, error = conditionMessage)
  stopifnot(grepl(pattern, message, fixed = TRUE))
}
fronts <- function(plot) {
  e <- export_ggchord_layout(plot, original_data = TRUE)
  d <- e$radial_endpoints
  f <- get_chord_layout(plot)$radial_frames
  ref <- get_chord_layout(plot)$sequence_reference
  # Every true interval front has inclusive width, even a one-base HSP.
  for (side in c("q", "s")) {
    ids <- d[[paste0(side, "accver")]]
    delta <- cbind(d[[paste0(side, "x_end")]] - d[[paste0(side, "x_start")]],
      d[[paste0(side, "y_end")]] - d[[paste0(side, "y_start")]])
    width <- sqrt(rowSums(delta^2))
    source_width <- abs(d[[paste0(side, "end")]] - d[[paste0(side, "start")]]) + 1
    stopifnot(isTRUE(all.equal(width, source_width * ggchord:::ggchord_collinear_scale(ids, ref), check.attributes = FALSE)))
  }
  e
}
roundtrip <- function(p) {
  l <- get_chord_layout(p)
  ids <- if (is.null(l$display_instances)) l$seqs else l$display_instances$instance_id
  pos <- unname(l$sequence_reference$lens[l$seqs])
  m <- map_ggchord_position(p, ids, pos)
  back <- map_ggchord_position(p, ids, m[c("x", "y")], inverse = TRUE)
  stopifnot(isTRUE(all.equal(back$position, pos, check.attributes = FALSE)))
  invisible(ggplot_gtable(ggplot_build(p)))
}

# Real sources: declared PAF lengths, complete nucleotide and protein FASTA.
paf <- read_paf("examples/format-inputs/emales-pair.paf")
s <- unique(rbind(data.frame(accver = paf$qaccver, length = paf$qlen),
  data.frame(accver = paf$saccver, length = paf$slen)))
p <- ggchord(s, paf, validate = "none") + geom_seq() + geom_link_ribbon() +
  geom_seq_label() + coord_radial_synteny()
stopifnot(nrow(fronts(p)$radial_endpoints) == nrow(paf))
roundtrip(p)
message("Real two-spoke PAF: ", nrow(paf), " source intervals retained once")
s <- read_fasta_lengths(files = list.files("examples/fasta", full.names = TRUE))
h <- read_blast(files = list.files("examples/blastn-standard", pattern = "tsv$", full.names = TRUE))
p <- ggchord(s, h, validate = "none") + geom_seq() + geom_link_ribbon() +
  geom_seq_label() + coord_radial_synteny(start_angle = 70)
e <- fronts(p)
stopifnot(sum(e$connection_report$retained) == nrow(e$radial_endpoints), identical(p$ggchord$sources$relationship, h))
roundtrip(p)
message("Real four-spoke nucleotide BLAST: ", nrow(e$radial_endpoints), "/", nrow(h), " HSPs selected by cyclic adjacency")
s <- read_fasta_lengths("examples/protein/four-proteins.faa"); s$sequence_unit <- "aa"
h <- read_blast("examples/protein/all-vs-all-blastp.tsv"); h$sequence_unit <- "aa"
p <- ggchord(s, h, validate = "none") + geom_seq() + geom_link_ribbon() +
  geom_seq_label() + coord_radial_synteny(clockwise = FALSE)
e <- fronts(p); roundtrip(p)
stopifnot(e$metadata$sequence_unit == "aa", all(is.na(e$sequence_layout$units_per_bp)))
message("Real four-spoke protein BLAST: ", nrow(e$radial_endpoints), "/", nrow(h), " HSPs; aa retained")

# Explicit synthetic topology/direction/boundary fixtures, not scientific data.
s <- data.frame(accver = c("A", "B", "C", "D", "E"), length = c(1000, 600, 400, 800, 50))
h <- data.frame(source_hsp_id = paste0("h", 1:5), qaccver = c("A", "B", "C", "D", "A"),
  saccver = c("B", "C", "D", "A", "C"), qstart = c(1, 20, 40, 200, 10),
  qend = c(300, 200, 40, 500, 30), sstart = c(10, 300, 60, 500, 80),
  send = c(200, 100, 60, 900, 90), pident = c(80, 90, 95, 100, 90))
base <- ggchord(s[1:4, ], h, validate = "none") + geom_seq() + geom_link_ribbon() + geom_seq_label()
p <- base + coord_radial_synteny()
e <- fronts(p)
stopifnot(setequal(e$radial_endpoints$source_row, 1:4),
  identical(e$connection_report$reason[5], "non_adjacent"))
# No D-A HSP => no fake closure; N=2 never duplicates a source row.
open <- ggchord(s[1:4, ], h[1:3, ], validate = "none") + geom_seq() + geom_link_ribbon() + coord_radial_synteny()
stopifnot(nrow(fronts(open)$radial_endpoints) == 3L)
two <- ggchord(s[1:2, ], h[1, ], validate = "none") + geom_seq() + geom_link_ribbon() + coord_radial_synteny()
stopifnot(nrow(fronts(two)$radial_endpoints) == 1L)
for (n in c(3, 5)) {
  keep <- h$qaccver %in% s$accver[seq_len(n)] & h$saccver %in% s$accver[seq_len(n)]
  roundtrip(ggchord(s[seq_len(n), ], h[keep, ], validate = "none") + geom_seq() +
    geom_link_ribbon() + coord_radial_synteny(start_angle = 25, clockwise = FALSE))
}
# No geometry or hidden routing creates extra HSP observations.
all_pairs <- ggchord(s[1:4, ], h, validate = "none") + geom_seq(seq_orientation = c(A=-1,B=1,C=-1,D=1), length_mode = "equal") +
  geom_link_ribbon(adjacent_only = FALSE) + coord_radial_synteny()
stopifnot(nrow(fronts(all_pairs)$radial_endpoints) == 5L)
roundtrip(all_pairs)
# Full endpoint widths and reverse twist survive, default paths miss the hole.
g <- export_ggchord_layout(p, include = "ribbon")$ribbon
stopifnot(all(sqrt(g$x^2 + g$y^2) >= .6 - 1e-7))
# Shared feature, label, track and source-space statistics through repeated instances.
f <- data.frame(accver = "A", start = c(100, 400), end = c(300, 500), strand = c("+", "-"), label = c("first", "second"))
events <- data.frame(accver = "A", pos = c(100, 100, 400))
rep <- ggchord(s[1:2, ], h[1, ], feature_data = f, validate = "none") +
  geom_seq(seq_order = c(A1 = "A", B = "B", A2 = "A")) +
  geom_link_ribbon(between = c("A1", "B")) + geom_feature(position = position_strand()) +
  geom_feature_label(label_orientation = "horizontal") +
  geom_track(data = events, geom = "bar", track = position_track("count", width = .1)) +
  coord_radial_synteny()
er <- fronts(rep); roundtrip(rep)
stopifnot(nrow(er$instance_routing) == 1L, er$instance_routing$q_instance == "A1")
stat <- get_chord_layout(rep)$track_statistics
stopifnot(length(stat) == 1L, all(tapply(stat[[1]]$y, stat[[1]]$accver, sum) == 3),
  length(unique(unlist(stat[[1]]$source_rows))) == 3L)
# Focused origins and boundary rendering remain distinct from one-based positions.
focused <- focus_ggchord_data(s[1:2, ], h[1, ], loci = data.frame(accver = c("A", "B"), start = c(20, 10), end = c(900, 500)))
fp <- ggchord(focused$seq_data, focused$ribbon_data, validate = "none") + geom_seq(seq_orientation = -1) + geom_link_ribbon() + coord_radial_synteny()
m <- map_ggchord_position(fp, "A", 100, source = TRUE)
stopifnot(abs(map_ggchord_position(fp, "A", m[c("x", "y")], source = TRUE, inverse = TRUE)$position - 100) < 1e-7)
fail(ggplot_build(ggchord(s[1, ]) + geom_seq() + coord_radial_synteny()), "at least two")
fail(ggplot_build(ggchord(s[1:2, ], h[1, ], validate = "none") +
  geom_seq(seq_order = c(A1="A", B="B", A2="A")) + geom_link_ribbon() + coord_radial_synteny()), "ambiguous")
fail(base + coord_radial_synteny() + geom_link_ribbon(link_branch="query"), "shared branches")
# Stable IDs require authorization even without repeated sequence instances.
fail(get_chord_layout(ggchord(s[1:2, ], h[1, ], validate="none") + geom_seq() + geom_link_ribbon() +
  geom_link_ribbon() + coord_radial_synteny()), "allow_repeated")
# A no-HSP plot still builds with a physical band; no empty ribbon scale.
blank <- ggchord(s[1:2, ]) + geom_seq() + geom_track(data=events, geom="bar",
  track=position_track("empty-relations", width=grid::unit(3,"mm"))) +
  coord_radial_synteny(start_angle=35, ratio=1.4)
invisible(ggplot_gtable(ggplot_build(blank)))
fail(base + coord_radial_synteny() + geom_seq_scale(), "shared straight scale bars")
# Ribbon aggregation remains source-local before radial reconstruction.
aggregated <- ggchord(s[1:2, ], h[1, ], validate="none") + geom_seq() +
  stat_ribbon_density(bins=5, weight="count") + coord_radial_synteny()
roundtrip(aggregated)
message("Synthetic cyclic topology, full fronts, repeated instances, Track statistics, focus and errors passed")
