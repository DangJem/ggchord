# Shared website course snippets; original figure recipes remain unchanged.
# ---- coursecomparison ----
comparison <- ggchord(radial_seqs, radial_hsp, feature_data = radial_features) +
  geom_seq() +
  geom_seq_label() +
  geom_link_ribbon(adjacent_only = FALSE) +
  geom_feature(shape = "arrow", width = 0.025)
p_coursechord <- comparison + coord_chord()
p_courselinear <- comparison + coord_collinear()
p_courseradial <- comparison + coord_radial_synteny()
p_coursesankey <- comparison + coord_sankey()
# ---- coursecomparisonreport ----
comparison_plots <- list(
  Chord = p_coursechord, Collinear = p_courselinear,
  Radial = p_courseradial, Sankey = p_coursesankey
)
comparison_reports <- lapply(comparison_plots, function(plot) {
  export_ggchord_layout(plot, original_data = TRUE)
})
names(comparison_reports)
# ---- coursesankeygap ----
p_coursesankeygap <- p_foundationsankey + coord_sankey(column_gap = 0.8)
p_coursesankeygap
# ---- courseradialangle ----
p_courseradialangle <- p_courseradial +
  coord_radial_synteny(start_angle = 45)
p_courseradialangle
# ---- coursefusion ----
course_center <- ggchord(fusion_seqs |> slice(1)) +
  geom_seq(show.legend = FALSE) +
  coord_circular()
course_detail <- ggchord(fusion_seqs) +
  geom_seq(show.legend = FALSE) +
  geom_seq_label() +
  coord_collinear(seq_axis = "vertical")
p_coursefusion <- ggchord_fusion(course_center) +
  fusion_panel(course_detail, id = "detail", side = "right") +
  fusion_ribbon(
    fusion_links |> slice(1), qpanel = "center", spanel = "detail",
    fill = "#269B82", alpha = 0.35
  )
p_coursefusion
# ---- coursecircleinput ----
library(dplyr)
circle_seqs <- plasmid_example_pUC19 |> select(accver, length)
circle_origins <- tibble(
  accver = "pUC19", start = 867, end = 1455,
  strand = "-", type = "rep_origin", anno = "ori"
)
p_coursecircle <- ggchord(circle_seqs) +
  geom_seq() + geom_seq_center_label() + coord_circular()
# ---- coursecirclefeature ----
p_coursecirclefeature <- p_coursecircle +
  geom_feature_plasmid(data = circle_origins) +
  geom_feature_label(data = circle_origins,
    label_side = "inside", label_radial_offset = 0.12, size = 4)
p_coursecirclefeature
# ---- courseblockb ----
library(ggchord)
library(dplyr)
# ---- courseblockc ----
ggsave("phage-comparison.png", p_real_chord, width = 10, height = 8, units = "in", dpi = 144)
# ---- courseblockd ----
ggsave("four-genome-chord.png", p_show_chord, width = 12, height = 10, units = "in", dpi = 144)
# ---- courseblocke ----
ggsave("ampR-neighbourhood.png", p_homologs, width = 10, height = 5, units = "in", dpi = 144)
# ---- courseblockf ----
ggsave("four-genome-collinear.png", p_show_collinear, width = 12, height = 8, units = "in", dpi = 144)
# ---- courseblockg ----
library(ggchord)
library(dplyr)
# ---- courseblockh ----
sequence <- plasmid_example_pSB1C3
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
# ---- courseblocki ----
sequence <- plasmid_example_pUC19
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
# ---- courseblockj ----
sequence <- plasmid_example_pBR322
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
# ---- courseblockk ----
sequence <- plasmid_example_pBluescript_II_SK_plus
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
# ---- courseblockl ----
sequence <- plasmid_example_pET_28a_plus
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
# ---- courseblockm ----
sequence <- plasmid_example_pETDuet_1
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
# ---- courseblockn ----
sequence <- plasmid_example_pDONR221
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
# ---- courseblocko ----
sequence <- plasmid_example_pcDNA3_1_plus
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
# ---- courseblockp ----
sequence <- plasmid_example_pTRE_Tight_BI
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
# ---- courseblockq ----
sequence <- plasmid_example_pTRIPZ
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
# ---- courseblockr ----
sequence <- plasmid_example_pCAMBIA1300
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
# ---- courseblocks ----
sequence <- plasmid_example_pEarleyGate_201
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
# ---- courseblockt ----
sequence <- plasmid_example_pSpCas9_BB_2A_GFP_PX458
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
# ---- courseblocku ----
ggsave("pUC19-map.png", p_plasmid, width = 10, height = 9, units = "in", dpi = 144)
# ---- courseblockv ----
ggsave("pUC19-tracks.png", p_multitrack, width = 9, height = 8, units = "in", dpi = 144)
# ---- courseblockw ----
library(ggchord)
library(dplyr)
seqs <- seq_data_example |>
  select(accver, length)
plot <- ggchord(seqs, ribbon_data_example) +
  geom_seq() +
  geom_link_ribbon() +
  adapt_ggchord(
    mode = "adaptive",
    reference_size = c(8, 6),
    text_scale_range = c(0.85, 1.35),
    min_text_size = 2.5
  )
report <- diagnose_ggchord(plot, width = 8, height = 6, device = "png")
report$guide_layout$guides
report$adaptive_layout$styles
# ---- courseblockx ----
data(plasmid_example_pUC19)
library(dplyr)
source_sequence <- plasmid_example_pUC19
bin_start <- seq.int(1L, source_sequence$length, by = 100L)
track_bins <- tibble(
  accver = source_sequence$accver, start = bin_start,
  end = pmin(bin_start + 99L, source_sequence$length)
) |>
  rowwise() |>
  mutate(gc = mean(strsplit(
    substr(source_sequence$sequence, start, end), "", fixed = TRUE
  )[[1]] %in% c("G", "C"))) |>
  ungroup()
themed_track <- ggchord(source_sequence |> select(accver, length)) +
  geom_seq() +
  geom_track(aes(x = start, y = gc), data = track_bins, geom = "point",
    track = position_track("GC", width = unit(5, "mm")), colour = "#59B59C") +
  scale_track_value_continuous("GC", limits = c(0, 1), unit = "fraction",
    guide = ggplot2::guide_axis(theme = theme(
      axis.text = element_text(colour = "#E5A5C0", size = 10),
      axis.title = element_text(face = "bold")
    ))) +
  coord_circular() +
  theme_ggchord_dark() +
  theme(
    ggchord.axis.text.orientation = "horizontal",
    ggchord.track.axis.text = element_text(size = 9),
    ggchord.track.axis.title.offset = unit(6, "mm")
  )
# ---- courseblocky ----
library(ggchord)
library(dplyr)
# ---- courseblockz ----
geom_feature_label_repel(label_segment_layer = "above")
# ---- courseblockaa ----
features <- find_common_features(plasmid_example_pTRE_Tight_BI)
feature_position <- position_feature_stack(base_position = position_plasmid())
ggchord(plasmid_example_pTRE_Tight_BI[c("accver", "length")], feature_data = features) +
  geom_seq(seq_style = "double") +
  geom_feature_plasmid(position = feature_position) +
  geom_feature_label_repel(
    position = feature_position,
    label_repeat = "bracket",
    label_fit = "shrink",
    label_min_size = 2.5
  ) +
  coord_circular() +
  theme_ggchord_plasmid()
# ---- courseblockab ----
geom_feature_label_repel(
  label_segment_colour = "grey50",
  label_segment_linetype = "dashed"
)
# ---- courseblockac ----
library(ggchord)
# ---- courseblockad ----
p_chord_feature + coord_collinear()
# ---- courseblockae ----
get_chord_layout(p_foundationsankey)$display_instances
# ---- courseblockaf ----
library(ggchord)
library(dplyr)
# ---- courseblockag ----
library(ggchord)
library(dplyr)
# ---- courseblockah ----
library(ggchord)
library(dplyr)
# ---- courseblockai ----
library(dplyr)
plasmid_example_pUC19 |> select(accver, length)
find_common_features(plasmid_example_pUC19) |> select(accver, start, end, strand, type, anno) |> slice_head(n = 5)
# ---- courseblockaj ----
library(ggchord)
library(dplyr)
# ---- courseblockak ----
library(ggchord)
library(dplyr)
# ---- courseblockal ----
ggsave("pETDuet-1-map.png", p_show_plasmid, width = 16.5, height = 8, units = "in", dpi = 144)
# ---- courseblockam ----
library(ggchord)
library(dplyr)
# ---- courseblockan ----
library(ggchord)
library(dplyr)
# ---- courseblockao ----
library(ggchord)
library(dplyr)
# ---- courseblockap ----
data(ribbon_data_example)
paired_rows <- ribbon_data_example |>
  transmute(
    chrom1 = qaccver, start1 = pmin(qstart, qend) - 1, end1 = pmax(qstart, qend),
    chrom2 = saccver, start2 = pmin(sstart, send) - 1, end2 = pmax(sstart, send)
  )
paired_file <- tempfile(fileext = ".bedpe")
write.table(paired_rows, paired_file,
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
paired_links <- read_bedpe(paired_file)
unlink(paired_file)
paired_links |>
  select(qaccver, saccver, qstart, qend, sstart, send, q_mapped, s_mapped)
# ---- courseblockaq ----
library(ggchord)
library(dplyr)
# ---- courseblockar ----
character_source <- read_fasta_lengths(
  "examples/plasmid/pUC19.fna",
  sequence_unit = "bp", keep_sequence = TRUE, source_file = TRUE
)
character_locus <- focus_ggchord_data(
  character_source,
  loci = data.frame(accver = character_source$accver, start = 1001, end = 1060)
)
character_locus$seq_data |>
  select(accver, length, sequence, .source_origin)
# ---- courseblockas ----
aligned_dna <- read_aligned_fasta(
  "examples/alignment/opuntia-aligned.fasta", "DNA", source_file = TRUE
)
selected_rows <- unique(aligned_dna$accver)[c(1, 3, 7)]
local_dna <- aligned_dna |>
  filter(accver %in% selected_rows, between(alignment_column, 321, 380))

map_alignment_position(local_dna, selected_rows[1], c(342, 343, 344))
map_alignment_position(local_dna, selected_rows[1], 343, inverse = TRUE)

p_characters <- ggplot2::ggplot(
  local_dna, aes(alignment_column, accver, label = symbol)
) +
  geom_residue(aes(fill = symbol), size = 3, overflow = "error") +
  ggplot2::scale_fill_manual(
    values = c(
      A = "#B8D8BA", C = "#A6CBE3", G = "#F2D28C", T = "#EAA8A8",
      N = "#D8D8D8", "-" = "#F5F5F5", "." = "#F5F5F5"
    ),
    guide = "none"
  ) +
  ggplot2::labs(x = "Alignment column", y = NULL) +
  ggplot2::theme_minimal()
p_characters
# ---- courseblockat ----
aligned_protein <- read_aligned_fasta(
  "examples/alignment/protein-aligned.fasta", "AA"
)
protein_rows <- unique(aligned_protein$accver)[2:4]
local_protein <- aligned_protein |>
  filter(accver %in% protein_rows, between(alignment_column, 221, 260))

ggplot2::ggplot(local_protein, aes(alignment_column, accver, label = symbol)) +
  geom_residue(aes(fill = symbol), size = 3, overflow = "tile") +
  ggplot2::scale_fill_hue(l = 85, c = 35, guide = "none") +
  ggplot2::labs(x = "Alignment column", y = NULL) +
  ggplot2::theme_minimal()
# ---- courseblockau ----
original_strings <- read_fasta_lengths(
  "examples/alignment/opuntia-aligned.fasta", keep_sequence = TRUE
)
fragment_sources <- original_strings |>
  transmute(
    accver, source_accver = accver, start = 1, strand = "+",
    sequence = gsub("[-.]", "", sequence)
  )
checked_cells <- as_ggchord_alignment(
  original_strings, "DNA", source_data = fragment_sources
)
# ---- courseblockav ----
library(dplyr)
# Requires the source checkout, dplyr and optional Biostrings.
genome <- read_fasta_lengths(
  "examples/format-inputs/emales-BVI_023A.fna", keep_sequence = TRUE
)
coding <- read_gff3("examples/format-inputs/emales-BVI_023A.gff3", unstranded = "keep") |>
  filter(feature_id == "BVI_023A_0009")
cds <- build_ggchord_cds(coding, genetic_code = 11)
observed <- translate_ggchord_cds(cds, genome, initiation = "start")
expected <- read_fasta_lengths(
  "examples/format-inputs/protein-pair.faa", keep_sequence = TRUE
) |>
  filter(accver == "BVI_023A_0009")
stopifnot(sub("[*]$", "", observed$models$protein_sequence) == expected$sequence)
observed$residues |>
  select(alignment_column, symbol, codon, source_positions) |>
  head()
# ---- courseblockaw ----
library(ggchord)
library(dplyr)
# ---- courseblockax ----
signal_file <- system.file("extdata", "pUC19-gc-100.bedGraph", package = "ggchord")
signal <- read_bedgraph(signal_file, source_file = TRUE)
data(plasmid_example_pUC19)
gc_band <- position_track("GC", side = "inside", width = unit(7, "mm"))
p_signal <- ggchord(plasmid_example_pUC19) +
  geom_seq(show.legend = FALSE) +
  geom_track(
    aes(xmin = start - 1, xmax = end, ymin = 0, ymax = value, fill = value),
    data = signal, geom = "rect", track = gc_band
  ) +
  scale_track_value_continuous("GC", limits = c(0, 100), unit = "GC (%)") +
  scale_track_fill_gradientn(
    "GC", colours = c("#F1E51D", "#31849B"), limits = c(0, 100), name = "GC (%)"
  ) +
  coord_circular()
p_signal

if (.Platform$OS.type != "windows" && requireNamespace("rtracklayer", quietly = TRUE)) {
  indexed_file <- system.file("extdata", "pUC19-gc-100.bw", package = "ggchord")
  selected_signal <- read_bigwig(
    indexed_file,
    which = data.frame(accver = "pUC19", start = c(50, 100), end = c(150, 250)),
    source_file = TRUE
  )
  selected_signal |> select(accver, start, end, value)
}
# ---- courseblockay ----
if (requireNamespace("GenomicAlignments", quietly = TRUE) &&
    requireNamespace("Rsamtools", quietly = TRUE)) {
  bam_file <- system.file("extdata", "sm_treated1.bam", package = "GenomicAlignments")
  read_filter <- Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(
      isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
      isNotPassingQualityControls = FALSE, isDuplicate = FALSE
    ), mapqFilter = 0L
  )
  aligned <- GenomicAlignments::readGAlignments(bam_file, param = read_filter)
  read_coverage <- GenomicAlignments::coverage(aligned, drop.D.ranges = TRUE)
  coverage_track <- bind_rows(lapply(names(read_coverage), function(id) {
    runs <- read_coverage[[id]]
    ends <- cumsum(S4Vectors::runLength(runs))
    data.frame(
      accver = id, start = c(1, head(ends, -1L) + 1), end = ends,
      value = S4Vectors::runValue(runs)
    )
  }))
  coverage_track |> head()
}
# ---- courseblockaz ----
library(ggchord)
library(dplyr)
# ---- courseblockba ----
if (requireNamespace("VariantAnnotation", quietly = TRUE)) {
  variant_file <- system.file("extdata", "chr22.vcf.gz", package = "VariantAnnotation")
  imported_variants <- read_vcf(variant_file, source_file = TRUE)
  imported_variants$features$variant |>
    select(accver, start, end, feature_id, REF, ALT) |>
    head()
}

if (requireNamespace("StructuralVariantAnnotation", quietly = TRUE)) {
  structural_file <- system.file("extdata", "gridss.vcf", package = "StructuralVariantAnnotation")
  structural <- read_vcf(structural_file, structural = TRUE)
  structural$links |>
    select(qaccver, saccver, qstart, qend, sstart, send,
      q_breakend_orientation, s_breakend_orientation)
  structural$report$records[[1]]$warnings
}
# ---- courseblockbb ----
alignment_python <- Sys.getenv("ALIGNMENT_PYTHON", "python3")
alignment_available <- nzchar(Sys.which(alignment_python)) &&
  identical(suppressWarnings(system2(
    alignment_python,
    c("-c", shQuote("from Bio import Align; assert hasattr(Align, 'parse')")),
    stdout = FALSE, stderr = FALSE
  )), 0L)
if (alignment_available) {
  maf <- read_maf(
    "examples/genome-alignments/ucsc_mm9_chr10.maf",
    alphabet = "DNA", python = alignment_python, source_file = TRUE
  )
  local_maf <- maf$alignment |>
    filter(.source_input == 1, source_block == 1, between(alignment_column, 61, 100))
  ggplot2::ggplot(local_maf, aes(alignment_column, accver, label = symbol)) +
    geom_residue(aes(fill = symbol), size = 3, overflow = "tile") +
    ggplot2::scale_fill_manual(
      values = c(
        A = "#B8D8BA", C = "#A6CBE3", G = "#F2D28C", T = "#EAA8A8",
        N = "#D8D8D8", "-" = "#F5F5F5"
      ), guide = "none"
    ) +
    ggplot2::labs(x = "MAF alignment column", y = NULL) +
    ggplot2::theme_minimal()

  chains <- read_chain(
    "examples/genome-alignments/dna_rna.chain", python = alignment_python
  )
  observed_chain <- chains$ribbon_data |>
    filter(source_block == 1)
  observed_chain |>
    select(qaccver, saccver, qstart, qend, sstart, send)
}
# ---- courseblockbc ----
delta_decoder <- Sys.getenv("DELTA_DECODER", "show-coords")
if (nzchar(Sys.which(delta_decoder))) {
  delta <- read_delta(
    "examples/genome-alignments/promer.delta",
    executable = delta_decoder, source_file = TRUE
  )
  observed_intervals <- delta$ribbon_data |>
    filter(qaccver == "3210101")
  displayed_ids <- unique(c(observed_intervals$qaccver, observed_intervals$saccver))
  displayed_sources <- delta$seq_data |>
    filter(accver %in% displayed_ids)

  ggchord(displayed_sources) +
    geom_seq(show.legend = FALSE) +
    geom_seq_label() +
    geom_link_ribbon(data = observed_intervals) +
    coord_collinear() +
    labs(
      title = "Existing translated alignment intervals",
      subtitle = "Drosophila source coordinates in bp; native identity in AA"
    )
}
# ---- courseblockbd ----
if (nzchar(Sys.which(delta_decoder))) {
  plant_delta <- read_delta(
    "examples/genome-alignments/arabidopsis-first-pair.delta",
    executable = delta_decoder, source_file = TRUE
  )
  plant_delta$ribbon_data |>
    select(qaccver, qstart, qend, saccver, sstart, send, pident) |>
    head()
}
# ---- courseblockbe ----
signature_file <- system.file("extdata", "P00517-signatures.tsv", package = "ggchord")
signatures <- read_signature_matches(signature_file, source_file = TRUE)
signatures$features$signature |>
  select(accver, start, end, signature_id, analysis, source_score)

if (requireNamespace("jsonlite", quietly = TRUE)) {
  cached_file <- system.file("extdata", "P51587-pfam.json", package = "ggchord")
  cached <- read_signature_matches(cached_file, format = "entry_json")
  ggchord(cached$seq_data) +
    geom_seq(show.legend = FALSE) +
    geom_seq_label() +
    geom_feature(
      aes(feature_fill = signature_id), data = cached$features$signature,
      shape = "block", position = position_feature_stack()
    ) +
    coord_collinear()
}
# ---- courseblockbf ----
pair_file <- system.file("extdata", "plant-feature-pairs.homology", package = "ggchord")
group_file <- system.file("extdata", "plant-feature-groups.txt", package = "ggchord")
position_file <- system.file("extdata", "plant-pair-features.gff", package = "ggchord")

# This four-column tool table is not GFF3. Explicitly declare source anchors.
anchors <- read.table(
  position_file, header = FALSE, na.strings = character(),
  col.names = c("accver", "feature_id", "start", "end"),
  colClasses = c("character", "character", "numeric", "numeric")
)
relations <- read_feature_relations(
  pair_file, feature_data = anchors, source_file = TRUE
)
relations$ribbon_data |>
  select(qfeature_id, sfeature_id, qaccver, qstart, qend, source_score)

groups <- read_feature_relations(group_file, format = "named_groups")
groups$membership |>
  select(group_id, feature_id, .source_group, .source_member)
# ---- courseblockbg ----
library(ggchord)
library(dplyr)
# ---- courseblockbh ----
p_count_shortcut <- ggchord(foundation_seqs) +
  geom_seq() +
  geom_track_bar(
    aes(x = pos), data = foundation_events,
    track = position_track("query_starts")
  ) +
  coord_collinear()
# ---- courseblockbi ----
estimate_bins <- bins |> dplyr::mutate(pos = (start + end) / 2)
estimate_band <- position_track("gc_fit", width = unit(12, "mm"))
p_estimate <- ggchord(track_seq) +
  geom_seq() +
  geom_track(
    aes(x = pos, y = gc), data = estimate_bins, geom = "ribbon", stat = "smooth",
    method = "lm", formula = y ~ x, se = TRUE, position = "identity",
    fill = "#387D88", alpha = 0.2, track = estimate_band
  ) +
  geom_track(
    aes(x = pos, y = gc), data = estimate_bins, geom = "line", stat = "smooth",
    method = "lm", formula = y ~ x, se = FALSE, colour = "#285B65", track = estimate_band
  ) +
  scale_track_value_continuous("gc_fit", name = "GC fit", unit = "fraction") +
  coord_collinear()
p_estimate
# ---- courseblockbj ----
library(ggchord)
library(dplyr)
# ---- courseblockbk ----
library(ggchord)
library(dplyr)
# ---- courseblockbl ----
# Requires dplyr and the source checkout or tutorial input bundle.
source("examples/fusion-routing.R")
routed
# ---- courseblockbm ----
library(ggchord)
library(dplyr)
# ---- courseblockbn ----
data(seq_data_example)
sequences <- seq_data_example |>
  dplyr::select(accver, length)
base_panels <- ggchord(sequences) + geom_seq(show.legend = FALSE)
central <- base_panels + coord_chord()
mirrored <- base_panels + coord_sankey()
free_angle <- ggchord_fusion(central) +
  fusion_panel(mirrored, "mirror", angle = 31, distance = .8, size = c(.16, .22))
free_box <- ggchord_fusion(central, center_box = c(.02, .35, .2, .8)) +
  fusion_panel(mirrored, "mirror", box = c(.45, .98, .1, .9))
physical_box <- ggchord_fusion(central, center_box = unit(c(5, 85, 10, 100), "mm")) +
  fusion_panel(mirrored, "mirror", box = unit(c(100, 210, 10, 100), "mm"))
export_ggchord_fusion(physical_box, width = 12, height = 6)$panels
# ---- courseblockbo ----
library(ggchord)
guide_source <- ggchord(seq_data_example) + geom_seq(show.legend = FALSE)
guide_fusion <- ggchord_fusion(guide_source, layout = "adaptive") +
  fusion_panel(guide_source + coord_collinear(), "detail") +
  fusion_ribbon(ribbon_data_example, aes(fill = qaccver, alpha = pident),
    qpanel = "center", spanel = "detail") +
  ggplot2::guides(
    fill = ggplot2::guide_legend(position = "bottom", nrow = 1),
    alpha = ggplot2::guide_legend(position = "right")
  )
guide_fusion
# ---- courseblockbp ----
# Requires dplyr and the source checkout or tutorial input bundle.
source("examples/fusion-auto.R")
auto_fusion
# ---- courseblockbq ----
styled_guide_fusion <- guide_fusion +
  ggplot2::labs(
    title = "Source-aware composition",
    subtitle = "Independent genomic panels and native Guides",
    caption = "Source: bundled ggchord records",
    tag = "A"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.caption = element_text(size = 9),
    plot.margin = margin(5, 5, 5, 5, "mm"),
    plot.background = element_rect(fill = "#F3F5F7", colour = NA)
  )
styled_guide_fusion
# ---- courseblockbr ----
library(ggchord)
library(dplyr)
# ---- courseblockbs ----
plasmid <- plasmid_example_pETDuet_1
annotations <- find_common_features(plasmid)
selected_region <- data.frame(accver = plasmid$accver, start = 101, end = 1900)
local_source <- focus_ggchord_data(
  plasmid,
  feature_data = annotations |> filter(start <= end),
  loci = selected_region
)
whole_plasmid <- ggchord(plasmid, feature_data = annotations) +
  geom_seq(show.legend = FALSE) +
  geom_feature(position = position_plasmid(), show.legend = FALSE) +
  geom_seq_center_label() +
  coord_circular(rotation = -30)
local_labels <- local_source$feature_data |>
  filter(type %in% c("CDS", "rep_origin"), anno != "ATG")
local_detail <- ggchord(
  local_source$seq_data, feature_data = local_source$feature_data
) +
  geom_seq(seq_orientation = -1, show.legend = FALSE) +
  geom_feature(show.legend = FALSE) +
  geom_feature_label(
    data = local_labels, label_orientation = "horizontal",
    label_overlap = "allow", colour = "#202020", size = 3
  ) +
  coord_collinear(seq_axis = "vertical") +
  scale_seq_position_continuous(coordinates = "source")
focused_view <- magnify_ggchord(whole_plasmid, selected_region, local_detail)
focused_view
focus_geometry <- export_ggchord_fusion(focused_view, width = 12, height = 8)
focus_geometry$focus_regions
focus_geometry$leaders
# ---- courseblockbt ----
library(dplyr)
seqs <- seq_data_example |> select(accver, length)
features <- gene_data_example |> mutate(strand = source_strand, type = 'CDS')
hsp <- ribbon_data_example |> slice(1)
regions <- bind_rows(
  transmute(hsp, accver = qaccver, start = pmin(qstart, qend), end = pmax(qstart, qend)),
  transmute(hsp, accver = saccver, start = pmin(sstart, send), end = pmax(sstart, send))
)
local <- focus_ggchord_data(seqs, feature_data = features,
  ribbon_data = hsp, loci = regions)
source <- ggchord(seqs, feature_data = features, ribbon_data = hsp) +
  geom_seq(show.legend = FALSE) + geom_gene(show.legend = FALSE) +
  geom_link_ribbon(show.legend = FALSE)
detail <- ggchord(local$seq_data, feature_data = local$feature_data,
  ribbon_data = local$ribbon_data) + geom_seq(show.legend = FALSE) +
  geom_gene(show.legend = FALSE) + geom_link_ribbon(show.legend = FALSE) + coord_collinear()
paired_detail <- magnify_ggchord(
  source + coord_collinear(), regions, detail, side = "bottom"
)
paired_detail
export_ggchord_fusion(paired_detail, width = 12, height = 8)$focus_regions
# ---- courseblockbu ----
library(ggchord)
library(dplyr)
seqs <- seq_data_example |> select(accver, length)
hsp <- ribbon_data_example |> slice(1)
regions <- bind_rows(
  transmute(hsp, accver = qaccver, start = pmin(qstart, qend),
    end = pmax(qstart, qend), detail_panel = "query"),
  transmute(hsp, accver = saccver, start = pmin(sstart, send),
    end = pmax(sstart, send), detail_panel = "subject")
)
make_detail <- function(region) {
  local <- focus_ggchord_data(seqs, loci = region)
  ggchord(local$seq_data) +
    geom_seq(show.legend = FALSE) +
    coord_collinear(seq_axis = "vertical")
}
source <- ggchord(seqs, hsp) +
  geom_seq(show.legend = FALSE) +
  geom_link_ribbon(show.legend = FALSE)
details <- list(
  query = make_detail(regions[1, ]),
  subject = make_detail(regions[2, ])
)
multiple <- magnify_ggchord(
  source, regions, details,
  side = c(query = "left", subject = "right")
)
export_ggchord_fusion(multiple, width = 12, height = 8)$focus_regions
# ---- courseblockbv ----
library(dplyr)
cells <- read_aligned_fasta("examples/alignment/opuntia-aligned.fasta", "DNA")
rows <- unique(cells$accver)[c(1, 3)]
window <- cells |>
  filter(accver %in% rows, between(alignment_column, 321, 380))
sequences <- cells |>
  filter(accver %in% rows, !is_gap) |>
  group_by(accver) |>
  summarise(length = max(source_position), .groups = "drop")
regions <- window |>
  filter(!is_gap) |>
  group_by(source_accver) |>
  summarise(start = min(source_position), end = max(source_position), .groups = "drop") |>
  rename(accver = source_accver)
letters <- ggplot2::ggplot(window, aes(alignment_column, accver, label = symbol)) +
  geom_residue(aes(fill = symbol), size = 2.7, overflow = "error") +
  ggplot2::scale_fill_hue(l = 85, c = 40, guide = "none") +
  ggplot2::theme_minimal()
letters <- register_ggchord_alignment(letters, window)
macro <- ggchord(sequences) + geom_seq() + coord_collinear()
character_detail <- magnify_ggchord(
  macro, regions, letters,
  center_size = c(.22, .7), detail_size = c(.5, .8)
)
export_ggchord_fusion(character_detail, width = 16, height = 5)$focus_regions
# ggsave("character-detail.png", character_detail, width = 16, height = 5)
# ---- courseblockbw ----
library(ggchord)
library(dplyr)
# ---- courseblockbx ----
data(plasmid_example_pUC19)
source_genome <- plasmid_example_pUC19 |>
  dplyr::select(accver, length)
large_region <- data.frame(accver = "pUC19", start = 101, end = 1900)
small_region <- data.frame(accver = "pUC19", start = 301, end = 600)
large_window <- focus_ggchord_data(source_genome, loci = large_region)
small_window <- focus_ggchord_data(source_genome, loci = small_region)
large_plot <- ggchord(large_window$seq_data) +
  geom_seq(show.legend = FALSE) + coord_collinear()
small_plot <- ggchord(small_window$seq_data) +
  geom_seq(show.legend = FALSE) + coord_collinear()
inner <- magnify_ggchord(large_plot, small_region, small_plot)
whole <- ggchord(source_genome) + geom_seq() + coord_circular()
large_region$detail_panel <- "detail/center"
outer <- magnify_ggchord(
  whole, large_region, inner,
  center_box = c(.02, .22, .15, .85), detail_box = c(.27, .99, .1, .9)
)
report <- export_ggchord_fusion(outer, width = 14, height = 7)
report$focus_regions
report$focus_region_data_by_composition
# ---- courseblockby ----
# Requires dplyr and the source checkout or tutorial input bundle.
source("examples/origin-detail.R")
origin_detail
# ---- courseend ----
