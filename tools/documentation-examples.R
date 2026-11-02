# Shared executable documentation snippets, current v0.16.6 development sources.
# Prefer source-backed bundled data; retain metadata outside plotting inputs.
# ---- setup ----
library(ggchord)
library(dplyr)
data(seq_data_example)
data(ribbon_data_example)
data(gene_data_example)
seqs <- seq_data_example |>
  select(accver, length)
links <- ribbon_data_example |>
  select(qaccver, saccver, qstart, qend, sstart, send)
# The fixture preserves biological strands separately from its old display strands.
genes <- gene_data_example |>
  transmute(accver, start, end, strand = source_strand)
# ---- sequencepalette ----
# Optional appearance choices; none of these columns or scales is required input.
library(ggchord)
seq_palette <- setNames(
  c("#3867D6", "#E58A28", "#219B83", "#9B59B6"),
  c("MT108731.1", "MT118296.1", "OQ646790.1", "OR222515.1")
)
gene_palette <- c("+" = "#3B6FD8", "-" = "#E1812C", "." = "#8291A6")
# ---- featurepalette ----
library(ggchord)
feature_palette <- c(
  CDS = "#3973CB", promoter = "#8B5EC4", rep_origin = "#E8A42D",
  primer_bind = "#168F86", misc_feature = "#A0ADBD", protein_bind = "#CC6485"
)
# ---- docstyle ----
library(ggchord)
doc_theme <- theme_ggchord_publication(
  axis.text = element_text(size = 11),
  axis.text.check.overlap = TRUE
) + theme(legend.position = "bottom")
doc_axis <- function(coordinates = "relative") {
  scale_seq_position_continuous(
    coordinates = coordinates, n.breaks = 3, minor_breaks = NULL,
    labels = function(x) format(x, scientific = FALSE, big.mark = ",", trim = TRUE)
  )
}
# ---- chord ----
p_chord <- ggchord(seq_data = seqs, ribbon_data = links, gene_data = genes) +
  geom_link_ribbon(aes(ribbon_fill = qaccver), alpha = .7) +
  scale_ribbon_fill_manual(values = seq_palette, guide = "none") +
  geom_seq(show.legend = FALSE) +
  geom_gene(position = position_strand(.1), width = .04) +
  geom_seq_label(size = 4, seq_label_radius = .84, seq_label_hjust = .5) +
  coord_chord() +
  scale_seq_colour_manual(values = seq_palette) +
  scale_gene_fill_manual(values = gene_palette) +
  doc_theme
p_chord
# ---- plasmiddata ----
library(dplyr)
data(plasmid_example_pUC19)
plasmid_seq <- plasmid_example_pUC19 |>
  select(accver, length)
feature_source <- find_common_features(plasmid_example_pUC19)
features <- feature_source |>
  select(accver, start, end, strand, type, anno)
# ---- circular ----
circle_features <- features |>
  filter(type %in% c("CDS", "rep_origin"))
p_circular <- ggchord(seq_data = plasmid_seq) +
  geom_seq() +
  geom_feature(aes(feature_fill = type), data = circle_features, position = position_plasmid()) +
  geom_feature_label_repel(data = circle_features, position = position_plasmid(), size = 4) +
  geom_seq_center_label() +
  scale_feature_fill_manual(values = feature_palette, breaks = unique(circle_features$type)) +
  coord_circular() +
  theme_ggchord_publication(axis.text = element_text(size = 11))
p_circular
# ---- collinear ----
anchors <- genes |>
  group_by(accver) |>
  slice_max(end - start, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(accver, start, end)
p_collinear <- ggchord(seq_data = seqs, ribbon_data = links, gene_data = genes) +
  geom_link_ribbon(aes(ribbon_fill = qaccver), alpha = .7) +
  scale_ribbon_fill_manual(values = seq_palette, guide = "none") +
  geom_seq() +
  geom_gene(position = position_strand(.1)) +
  geom_seq_label(size = 4) +
  coord_collinear(anchor = anchors, align = "center", orientation = "auto") +
  doc_axis() +
  scale_seq_colour_manual(values = seq_palette) +
  scale_gene_fill_manual(values = gene_palette) +
  doc_theme
p_collinear
# ---- adjustments ----
identity_links <- ribbon_data_example |>
  select(qaccver, saccver, qstart, qend, sstart, send, pident)
p_mapped <- ggchord(seq_data = seqs, ribbon_data = identity_links) +
  geom_seq() +
  geom_link_ribbon(aes(ribbon_fill = pident), alpha = .6) +
  scale_ribbon_fill_gradientn(colours = c("#3D2185", "#C74678", "#F6B739"), name = "Identity (%)") +
  scale_seq_colour_manual(values = seq_palette) +
  coord_chord(rotation = 20) +
  labs(title = "Sequence comparison") +
  doc_theme
p_mapped
# ---- labels ----
label_data <- gene_data_example |>
  filter(anno != "hypothetical protein") |>
  transmute(accver, start, end, strand = source_strand, anno)
p_labels <- ggchord(seq_data = seqs, gene_data = genes) +
  geom_seq() +
  geom_gene() +
  geom_gene_label_repel(data = label_data, label_layout = "auto", size = 4, label_wrap = 22) +
  geom_seq_label(size = 4) +
  coord_collinear() +
  doc_axis() +
  scale_seq_colour_manual(values = seq_palette) +
  scale_gene_fill_manual(values = gene_palette) +
  doc_theme
p_labels
# ---- trackdata ----
library(dplyr)
data(plasmid_example_pUC19)
track_seq <- plasmid_example_pUC19 |>
  select(accver, length)
starts <- seq.int(1L, track_seq$length, by = 100L)
gc_fraction <- function(sequence) {
  bases <- strsplit(sequence, "")[[1]]
  mean(bases %in% c("G", "C"))
}
bins <- tibble(
  accver = track_seq$accver,
  start = starts,
  end = pmin(starts + 99L, track_seq$length)
) |>
  rowwise() |>
  mutate(gc = gc_fraction(substr(plasmid_example_pUC19$sequence, start, end))) |>
  ungroup()
head(bins)
# ---- tracks ----
p_tracks <- ggchord(seq_data = track_seq) +
  geom_seq() +
  geom_seq_center_label() +
  geom_track_bar(aes(track_value = gc, track_fill = gc),
    data = bins,
    position = position_track("gc", width = unit(10, "mm"))
  ) +
  scale_track_value_continuous("gc", limits = c(0, 1), unit = "GC fraction") +
  scale_track_fill_gradientn("gc", colours = c("#E9F1FE", "#65A1DC", "#224B9B")) +
  coord_circular() +
  doc_theme
p_tracks
# ---- recognition ----
data(plasmid_example_pUC19)
target <- setNames(plasmid_example_pUC19$sequence, plasmid_example_pUC19$accver)
primers <- find_primer_bindings(target, set = "universal")
sites <- find_restriction_sites(target, enzymes = c("EcoRI", "BamHI", "HindIII"), methylation = "none")
recognition_seq <- plasmid_example_pUC19 |>
  select(accver, length)
p_sites <- ggchord(seq_data = recognition_seq) +
  geom_seq() +
  geom_primer(data = primers, colour = "#168F86") +
  geom_restriction_site(data = sites, colour = "#29364B") +
  coord_circular()
p_sites
# ---- homologinput ----
library(dplyr)
data(plasmid_example_pUC19)
data(plasmid_example_pETDuet_1)
homolog_source <- bind_rows(plasmid_example_pUC19, plasmid_example_pETDuet_1)
homolog_seq <- homolog_source |>
  select(accver, length)
cluster_features <- homolog_source |>
  find_common_features() |>
  filter(type == "CDS", anno == "AmpR") |>
  select(accver, start, end, strand, anno) |>
  arrange(match(accver, homolog_source$accver)) |>
  mutate(feature_id = paste(accver, start, end, sep = ":"))
# Membership is computed from exact strand-corrected CDS DNA equality, not gene names.
reverse_complement <- function(sequence) {
  bases <- strsplit(sequence, "")[[1]]
  chartr("ACGT", "TGCA", paste(rev(bases), collapse = ""))
}
coding_sequences <- cluster_features |>
  left_join(select(homolog_source, accver, sequence), by = "accver") |>
  rowwise() |>
  mutate(
    coding_dna = substr(sequence, start, end),
    coding_dna = if (strand == "-") reverse_complement(coding_dna) else coding_dna
  ) |>
  ungroup()
cds_dna <- pull(coding_sequences, coding_dna)
# ---- homologmembership ----
membership <- tibble(
  feature_id = cluster_features$feature_id,
  cluster_id = paste0("exact_CDS_", match(cds_dna, unique(cds_dna)))
)
homologs <- link_ggchord_clusters(cluster_features, membership,
  pairs = tibble(qaccver = homolog_seq$accver[1], saccver = homolog_seq$accver[2])
)
loci <- find_ggchord_loci(homolog_seq, cluster_features,
  targets = cluster_features$feature_id[cluster_features$anno == "AmpR"][1], upstream = 100, downstream = 100
)
# ---- homologs ----
# Compare equivalent local windows instead of shrinking one gene on full plasmids.
homolog_windows <- cluster_features |>
  select(accver, start, end)
local_homologs <- focus_ggchord_data(homolog_seq, homologs$data, cluster_features,
  loci = homolog_windows, expand = 200
)
p_homologs <- ggchord(
  seq_data = local_homologs$seq_data,
  ribbon_data = local_homologs$ribbon_data, gene_data = local_homologs$gene_data
) +
  geom_link_ribbon(fill = "#8D7AC6", alpha = .65) +
  geom_seq(colour = "#253954", show.legend = FALSE) +
  geom_feature(
    data = local_homologs$gene_data, feature_fill = "#D79A24",
    position = position_strand(.10), width = .07
  ) +
  geom_feature_label_repel(data = local_homologs$gene_data, size = 4, label_fit = "shrink") +
  geom_seq_label(size = 4) +
  geom_seq_break() +
  coord_collinear() +
  doc_axis("source") +
  doc_theme
p_homologs
# ---- projection ----
# Run the exact-CDS comparison first; these are complete source-backed CDS.
coding <- cluster_features |>
  mutate(phase = 0)
cds <- build_ggchord_cds(coding, genetic_code = 11)
# This selects residues for a display example; it does not claim a detected domain.
domains <- tibble(feature_id = "selected_residues", parent_id = coding$feature_id[1], start = 2, end = 20)
genomic <- project_ggchord_features(domains, cds)
genomic
# ---- focus ----
focus_loci <- tibble(
  accver = seqs$accver,
  start = floor(seqs$length * .25), end = floor(seqs$length * .75)
)
focused <- focus_ggchord_data(seqs, links, genes, loci = focus_loci)
p_focus <- ggchord(
  seq_data = focused$seq_data, ribbon_data = focused$ribbon_data,
  gene_data = focused$gene_data
) +
  geom_seq() +
  geom_link_ribbon(fill = "#7199DF", alpha = .5) +
  geom_gene() +
  geom_seq_break() +
  geom_seq_label(size = 4) +
  scale_seq_colour_manual(values = seq_palette) +
  scale_gene_fill_manual(values = gene_palette) +
  doc_axis("source") +
  coord_collinear() +
  doc_theme
p_focus
# ---- segment ----
windows <- tibble(
  accver = rep(seqs$accver[1], 2),
  start = c(1, floor(seqs$length[1] * .65)), end = c(floor(seqs$length[1] * .4), seqs$length[1]),
  locus_id = c("left_window", "right_window"), row_id = seqs$accver[1], orientation = c(1, -1)
)
segmented <- segment_ggchord_data(slice_head(seqs, n = 1), windows = windows)
p_segments <- ggchord(seq_data = segmented$seq_data) +
  geom_seq() +
  coord_collinear(segments = segmented$segments, segment_gap = 2000)
p_segments
# ---- validation ----
validation <- validate_ggchord_data(seqs, links, genes)
cleaned <- clean_ggchord_data(seqs, links, genes)
prepared <- prepare_ggchord_data(seqs, links, features = list(genes = genes))
# ---- provenance ----
geometry <- export_ggchord_layout(p_collinear, include = "seq")
rows <- geometry$seq |>
  distinct(layer_id, component, source_row)
# ---- migration ----
p_migrated <- ggchord(plasmid_seq) +
  geom_seq() +
  geom_feature(data = features, shape = "block", width = .05) +
  geom_feature_label(data = features, label_rotation = 0) +
  coord_circular()
p_offsets <- ggchord(seqs, gene_data = genes) +
  geom_seq() +
  geom_gene(position = position_strand(c("+" = .1, "-" = -.2)))
# ---- save ----
ggsave("comparison.png", p_collinear, width = 8, height = 5, units = "in", dpi = 144)
# ---- preview ----
view_ggchord(p_collinear, width = 8, height = 5, units = "in")
diagnosis <- diagnose_ggchord(p_collinear, width = 8, height = 5, units = "in")
# ---- imports ----
seq_files <- read_fasta_lengths(files = "examples/fasta/*.fna", source_file = TRUE)
link_files <- read_blast(files = "examples/blastn-standard/*.tsv", source_file = TRUE)
feature_files <- read_gff3(files = "examples/gff3/*.gff3", unstranded = "keep", source_file = TRUE)
# ---- realimport ----
library(dplyr)
sequence_source <- read_fasta_lengths(files = "examples/fasta/*.fna", source_file = TRUE)
seq_real <- sequence_source |>
  select(accver, length)
annotation_source <- read_gff3(files = "examples/gff3/*.gff3", unstranded = "keep", source_file = TRUE)
feature_real <- annotation_source |>
  select(accver, start, end, strand, type, anno)
alignment_source <- read_blast(files = "examples/blastn-standard/*.tsv", source_file = TRUE)
link_real <- alignment_source |>
  select(qaccver, saccver, qstart, qend, sstart, send, pident, length)
real_validation <- validate_ggchord_data(seq_real, link_real)
real_validation
# ---- realprepare ----
real_prepared <- prepare_ggchord_data(seq_real, link_real,
  features = list(annotation = feature_real)
)
real_filtered <- filter_ggchord_ribbons(link_real, min_pident = 80, min_length = 500)
# Explicit display subset: the six longest CDS per sequence; raw annotations stay intact.
real_display_features <- feature_real |>
  group_by(accver) |>
  slice_max(end - start + 1, n = 6, with_ties = FALSE) |>
  ungroup()
real_labels <- real_display_features |>
  filter(anno != "hypothetical protein") |>
  group_by(accver) |>
  slice_head(n = 2) |>
  ungroup()
real_filtered$report
# ---- realchord ----
p_real_chord <- ggchord(seq_data = seq_real, ribbon_data = real_filtered$data, validate = "none") +
  geom_link_ribbon(aes(ribbon_fill = pident), alpha = .65) +
  geom_seq(show.legend = FALSE) +
  geom_feature(aes(feature_fill = strand),
    data = real_display_features,
    width = .045, position = position_strand(.08)
  ) +
  geom_feature_label_repel(
    data = real_labels, label_layout = "radial",
    label_wrap = 20, size = 4, position = position_strand(.08)
  ) +
  geom_seq_label(size = 4, seq_label_radius = .84, seq_label_hjust = .5) +
  coord_chord() +
  scale_seq_colour_manual(values = seq_palette) +
  scale_feature_fill_manual(values = gene_palette, breaks = c("+", "-")) +
  scale_ribbon_fill_gradientn(colours = c("#3D2185", "#C74678", "#F6B739"), name = "Identity (%)") +
  doc_theme
p_real_chord
# ---- reallinear ----
p_real_linear <- ggchord(seq_real, real_filtered$data, validate = "none") +
  geom_link_ribbon(fill = "#7199DF", alpha = .45) +
  geom_seq() +
  geom_seq_label(size = 4) +
  geom_feature(data = feature_real, feature_fill = "#3B6FD8", show.legend = FALSE, position = position_strand(.18)) +
  geom_seq_scale(length = 10000) +
  coord_collinear() +
  doc_axis() +
  scale_seq_colour_manual(values = seq_palette) +
  doc_theme
p_real_linear
# ---- plasmidinput ----
library(dplyr)
data(plasmid_example_pUC19)
plasmid_seq <- plasmid_example_pUC19 |>
  select(accver, length)
plasmid_feature_source <- find_common_features(plasmid_example_pUC19)
# Dedicated primer layers already show bindings; avoid duplicate primer annotations.
plasmid_features <- plasmid_feature_source |>
  filter(type != "primer_bind") |>
  select(accver, start, end, strand, type, anno)
plasmid_primers <- find_primer_bindings(plasmid_example_pUC19, set = "universal")
plasmid_sites <- find_restriction_sites(plasmid_example_pUC19)
display_sites <- filter_restriction_sites(plasmid_sites, enzymes = c("EcoRI", "BamHI", "HindIII"))
primer_names <- primer_catalog(set = "universal")
# ---- plasmid ----
plasmid_position <- position_feature_stack(spacing = .10, base_position = position_plasmid())
p_plasmid <- ggchord(seq_data = plasmid_seq) +
  geom_seq(seq_style = "double") +
  geom_feature_plasmid(aes(feature_fill = type), data = plasmid_features, position = plasmid_position) +
  geom_feature_label_repel(
    data = plasmid_features, position = plasmid_position,
    size = 4
  ) +
  geom_primer(data = plasmid_primers, colour = "#168F86") +
  geom_primer_label_repel(data = plasmid_primers, colour = "#167F80", max_overlaps = 0, size = 3) +
  geom_restriction_site(data = display_sites, colour = "#29364B", label_size = 3, leader = "radial") +
  geom_seq_center_label() +
  scale_feature_fill_manual(values = feature_palette) +
  scale_seq_position_continuous() +
  coord_circular(rotation = 90) +
  theme_ggchord_plasmid()
p_plasmid
# ---- cleaning ----
# Deliberately corrupt copies to demonstrate diagnostics; original rows remain intact.
problem_links <- ribbon_data_example |>
  slice_head(n = 2) |>
  select(qaccver, saccver, qstart, qend, sstart, send, pident) |>
  mutate(
    qaccver = replace(qaccver, 2, "absent"),
    pident = replace(pident, 1, 103)
  )
problem_report <- validate_ggchord_data(seqs, problem_links)
repaired <- clean_ggchord_data(seqs, problem_links, unknown_id = "drop", invalid_pident = "clip")
# ---- hierarchy ----
# Each source record represents a separate complete phage genome.
contigs <- seqs |>
  mutate(genome = accver)
genomes <- tibble(genome = seqs$accver)
hierarchy <- validate_ggchord_hierarchy(contigs, genomes, strict = TRUE)
# ---- adapters ----
renamed_sequences <- seqs |>
  rename(chrom = accver, size = length)
sequence_input <- as_ggchord_seq(renamed_sequences, mapping = c(accver = "chrom", length = "size"))
renamed_features <- genes |>
  rename(chrom = accver, left = start, right = end)
feature_input <- as_ggchord_feature(renamed_features,
  mapping = c(accver = "chrom", start = "left", end = "right")
)
link_input <- as_ggchord_link(links)
# ---- formats ----
# BED3 carries intervals only; the source annotations remain in features.
bed_path <- tempfile(fileext = ".bed")
bed_rows <- tibble(accver = features$accver, start0 = features$start - 1L, end = features$end)
write.table(bed_rows, bed_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
bed_features <- read_bed(bed_path)
unlink(bed_path)
# Unchanged PAF rows from the public EMALE tutorial, not a simulated alignment.
paf_source <- read_paf("examples/format-inputs/emales-pair.paf")
paf_selected <- paf_source |>
  filter(!is.na(mapq), mapq >= 30)
paf_links <- paf_selected |>
  select(qaccver, saccver, qstart, qend, sstart, send)
paf_seqs <- bind_rows(
  transmute(paf_source, accver = qaccver, length = qlen),
  transmute(paf_source, accver = saccver, length = slen)
) |>
  distinct()
p_paf <- ggchord(seq_data = paf_seqs, ribbon_data = paf_links) +
  geom_link_ribbon(fill = "#7199DF", alpha = .5) +
  geom_seq(colour = "#253954") +
  geom_seq_label(size = 4) +
  coord_collinear() +
  doc_axis() +
  doc_theme
p_paf
# ---- standardblast ----
blast_source <- read_blast("examples/blastn-standard/MT108731.1__MT118296.1.tsv")
blast_links <- blast_source |>
  select(qaccver, saccver, qstart, qend, sstart, send)
# A six-column table is sufficient for drawing; extra fields stay in blast_source.
custom_path <- tempfile(fileext = ".tsv")
write.table(blast_links, custom_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
custom_source <- read_blast(custom_path, format = "custom", col_names = names(blast_links))
custom_links <- custom_source |>
  select(all_of(names(blast_links)))
unlink(custom_path)
# ---- protein ----
protein_source <- read_blast("examples/format-inputs/emales-mavirus-blastp.tsv", source_file = TRUE)
emale_source <- read_gff3("examples/format-inputs/emales-BVI_023A.gff3", source_file = TRUE)
# The source declares complete Prodigal CDS and translation table 11.
emale_coding <- emale_source |>
  filter(grepl("partial=00", attributes)) |>
  select(feature_id, accver, start, end, strand, phase)
emale_cds <- build_ggchord_cds(emale_coding, genetic_code = 11)
# Keep strong HSPs whose query protein has a known complete CDS.
# Gaps are retained: projecting one query span does not assert an ungapped pair.
protein_kept <- protein_source |>
  filter(qaccver %in% emale_coding$feature_id, evalue <= 1e-10, length >= 50) |>
  arrange(desc(bitscore))
selected_hsp <- protein_kept |>
  slice_head(n = 1)
query_interval <- tibble(
  feature_id = "query_HSP", parent_id = selected_hsp$qaccver,
  start = selected_hsp$qstart, end = selected_hsp$qend
)
projected_query <- project_ggchord_features(query_interval, emale_cds)
# Subject projection additionally requires its own CDS map; it is not guessed from an accession.
# ---- genbank ----
# The original NCBI L09137.2 record identifies pUC19c; keep that identity.
genbank_result <- NULL
if (nzchar(Sys.which("python3")) && identical(suppressWarnings(system2("python3",
  c("-c", shQuote("from Bio import SeqIO")),
  stdout = FALSE, stderr = FALSE
)), 0L)) {
  genbank_result <- read_genbank("examples/format-inputs/pUC19c-L09137.2.gb")
}
# ---- optionalranges ----
range_inputs <- NULL
if (requireNamespace("GenomicRanges", quietly = TRUE)) {
  unknown <- features |>
    filter(strand == ".") |>
    slice_head(n = 1)
  ranges <- GenomicRanges::GRanges(
    seqnames = unknown$accver,
    ranges = IRanges::IRanges(unknown$start, unknown$end), strand = "*",
    seqlengths = setNames(plasmid_seq$length, plasmid_seq$accver)
  )
  range_inputs <- list(seq_data = as_ggchord_seq(ranges), features = as_ggchord_feature(ranges))
}
# ---- ribbons ----
unique_links <- deduplicate_ggchord_ribbons(real_filtered$data,
  by = "overlap",
  min_reciprocal_overlap = .9
)
merged_links <- merge_ggchord_ribbons(unique_links$data,
  max_gap = 100,
  min_pident_difference = 2
)
bundled_links <- bundle_ggchord_ribbons(link_real, seq_real, bins = 30)
p_density <- ggchord(seq_real, link_real, validate = "none") +
  geom_seq(show.legend = FALSE) +
  geom_seq_label(size = 4, seq_label_radius = .84, seq_label_hjust = .5) +
  stat_ribbon_density(aes(ribbon_alpha = after_stat(density)), fill = "#5379CE") +
  scale_ribbon_alpha_continuous(
    name = "Relative density", range = c(.1, .6),
    guide = guide_ggchord_legend(position = "bottom", nrow = 1)
  ) +
  scale_seq_colour_manual(values = seq_palette) +
  coord_chord() +
  doc_theme
p_bundle <- ggchord(seq_real, link_real, validate = "none") +
  geom_seq() +
  stat_ribbon_bundle(bins = 30, fill = "#5379CE") +
  coord_chord() +
  doc_theme
p_density
# ---- pointlinks ----
# Midpoints summarize existing HSPs for display; they are not new alignments.
points <- ribbon_data_example |>
  transmute(
    qaccver, saccver,
    qpos = floor((qstart + qend) / 2),
    spos = floor((sstart + send) / 2),
    pident
  )
p_points <- ggchord(seq_data = seqs) +
  geom_seq() +
  geom_seq_label(size = 4) +
  geom_link_line(aes(link_colour = pident),
    data = points,
    arrow = arrow(length = unit(2, "mm")), linewidth = .7
  ) +
  scale_link_colour_gradientn(
    colours = c("#3D2185", "#F6B739"),
    name = "Identity (%)", guide = guide_ggchord_colourbar(
      available_aes = "link_colour",
      position = "bottom", direction = "horizontal",
      theme = theme(legend.key.width = unit(35, "mm"), legend.key.height = unit(3, "mm"))
    )
  ) +
  scale_seq_colour_manual(values = seq_palette) +
  coord_chord() +
  doc_theme
p_points
# ---- featurestack ----
stack_features <- features |>
  select(accver, start, end, strand, type, anno)
stack_position <- position_feature_stack(spacing = .12)
feature_shapes <- setNames(ifelse(names(feature_palette) == "CDS", "arrow", "block"), names(feature_palette))
p_stack <- ggchord(seq_data = plasmid_seq) +
  geom_seq(show.legend = FALSE) +
  geom_feature(aes(feature_fill = type, feature_shape = type),
    data = stack_features,
    width = .08, position = stack_position
  ) +
  geom_feature_label_repel(data = stack_features, position = stack_position, size = 3.8, label_fit = "shrink", label_min_size = 2.8) +
  scale_feature_shape_manual(values = feature_shapes) +
  scale_feature_fill_manual(values = feature_palette) +
  coord_collinear() +
  doc_axis() +
  doc_theme
p_stack
# ---- labellayouts ----
# Compare layouts with at most the first and last informative label per sequence.
comparison_labels <- label_data |>
  group_by(accver) |>
  filter(row_number() %in% c(1L, n())) |>
  ungroup()
label_plots <- lapply(c("auto", "radial", "arc"), function(layout) {
  ggchord(seq_data = seqs, gene_data = genes) +
    geom_seq() +
    geom_gene() +
    geom_gene_label_repel(data = comparison_labels, label_layout = layout, size = 4, label_wrap = 22) +
    geom_seq_label(size = 4) +
    scale_seq_colour_manual(values = seq_palette) +
    scale_gene_fill_manual(values = gene_palette) +
    coord_chord(xlim = c(-1.8, 1.8), ylim = c(-1.6, 1.6)) +
    theme_ggchord_publication(axis = element_blank()) +
    labs(title = layout)
})
names(label_plots) <- c("auto", "radial", "arc")
p_manual_labels <- ggchord(seq_data = seqs, gene_data = genes) +
  geom_seq() +
  geom_gene() +
  geom_gene_label(data = label_data, label_rotation = 30) +
  coord_collinear()
label_plots$auto
# ---- decorations ----
regions <- tibble(
  accver = links$qaccver[1], start = min(links$qstart[1], links$qend[1]),
  end = max(links$qstart[1], links$qend[1]), category = "alignment interval"
)
p_decorated <- ggchord(seq_data = seqs, ribbon_data = links) +
  geom_seq_region(aes(region_fill = category), data = regions, colour = NA) +
  geom_link_ribbon(fill = "#7199DF", alpha = .4) +
  geom_seq() +
  geom_seq_label(size = 4) +
  scale_region_fill_manual(values = c("alignment interval" = "#E3D0A5")) +
  scale_seq_colour_manual(values = seq_palette) +
  geom_seq_scale(length = 10000) +
  coord_collinear() +
  doc_axis() +
  doc_theme
p_decorated
# ---- multitrack ----
line_bins <- bins |>
  transmute(
    accver,
    pos = floor((start + end) / 2),
    gc_mean = as.numeric(stats::filter(gc, rep(1 / 3, 3), sides = 2))
  ) |>
  filter(is.finite(gc_mean))
categories <- features |>
  select(accver, start, end, type)
p_multitrack <- ggchord(seq_data = track_seq) +
  geom_seq() +
  geom_seq_center_label() +
  geom_track_bar(aes(track_value = gc, track_fill = gc),
    data = bins,
    position = position_track("gc", width = unit(10, "mm"))
  ) +
  geom_track_line(aes(track_value = gc_mean),
    data = line_bins,
    position = position_track("gc", width = unit(10, "mm")), colour = "#454D65"
  ) +
  scale_track_value_continuous("gc", limits = c(0, 1), unit = "GC fraction") +
  scale_track_fill_gradientn("gc", colours = c("#E9F1FE", "#65A1DC", "#224B9B")) +
  geom_track_tile(aes(track_fill = type),
    data = categories,
    position = position_track("features", order = 2, width = unit(2, "mm"))
  ) +
  scale_track_fill_manual("features", values = feature_palette) +
  coord_circular() +
  doc_theme +
  theme(legend.box = "vertical", legend.box.just = "left")
focused_bins <- focus_ggchord_track(bins, track_seq,
  loci = tibble(accver = track_seq$accver, start = 500, end = 1500), type = "interval"
)
p_multitrack
# ---- linkprojection ----
# Exact CDS DNA equality supplies evidence for these selected ungapped residues.
matching <- homologs$data |>
  slice_head(n = 1)
paired_coding <- cluster_features |>
  mutate(phase = 0)
paired_cds <- build_ggchord_cds(paired_coding)
protein_link <- tibble(
  qfeature_id = matching$qfeature_id, sfeature_id = matching$sfeature_id,
  qstart = 1, qend = 20, sstart = 1, send = 20
)
projected_links <- project_ggchord_links(protein_link, paired_cds, ungapped = TRUE)
projected_links
# ---- mappingposition ----
source_position <- focus_loci$start[1] + 100
source_mark <- map_ggchord_position(p_focus, seqs$accver[1], source_position, source = TRUE)
p_marked <- p_focus +
  ggplot2::geom_point(
    data = source_mark, aes(x, y),
    inherit.aes = FALSE, colour = "#9B59B6", size = 3
  )
p_marked
# ---- styling ----
p_styled <- ggchord(seq_data = seqs, ribbon_data = links, gene_data = genes) +
  geom_seq() +
  geom_link_ribbon(aes(ribbon_fill = qaccver), alpha = .5) +
  geom_gene() +
  geom_seq_label(size = 4) +
  coord_collinear() +
  doc_axis() +
  scale_seq_colour_manual(values = seq_palette) +
  scale_gene_fill_manual(values = gene_palette) +
  scale_ribbon_fill_manual(values = seq_palette) +
  guides(gene_fill = guide_ggchord_legend(nrow = 1)) +
  doc_theme
theme_variants <- list(minimal = theme_ggchord_minimal(), dark = theme_ggchord_dark(), plasmid = theme_ggchord_plasmid())
p_styled
# ---- layoutresults ----
full_layout <- get_chord_layout(p_collinear)
layer_tables <- export_ggchord_layout(p_collinear, include = c("seq", "gene"), original_data = TRUE)
layout_choice <- optimize_ggchord_layout(seqs, links, optimize = c("order", "orientation"), weight = "count")
# ---- parameters ----
# Sequence IDs remain stable when the data is reordered; positional lists do not.
radii <- setNames(c(1.1, 1, 1, 1), seqs$accver)
widths <- c("+" = .045, "-" = .065)
p_parameters <- ggchord(seq_data = seqs, gene_data = genes) +
  geom_seq(seq_radius = radii) +
  geom_gene(
    width = widths,
    position = position_strand(c("+" = .08, "-" = -.1))
  ) +
  geom_seq_label(size = 4) +
  scale_seq_colour_manual(values = seq_palette) +
  scale_gene_fill_manual(values = gene_palette) +
  coord_chord() +
  doc_theme
p_parameters
# ---- firstdata ----
library(ggchord)
library(dplyr)
data(seq_data_example)
data(gene_data_example)
data(ribbon_data_example)
ids <- c("MT108731.1", "MT118296.1")
first_sequences <- seq_data_example |>
  filter(accver %in% ids) |>
  select(accver, length)
first_genes <- gene_data_example |>
  filter(accver %in% ids) |>
  transmute(accver, start, end, strand = source_strand)
first_links <- ribbon_data_example |>
  filter(qaccver %in% ids, saccver %in% ids) |>
  select(qaccver, saccver, qstart, qend, sstart, send)
# ---- firstplot ----
first_colours <- setNames(c("#3264D6", "#E18428"), ids)
p_first <- ggchord(
  seq_data = first_sequences, ribbon_data = first_links,
  gene_data = first_genes
) +
  geom_link_ribbon(fill = "#A4B3CC", alpha = .55) +
  geom_seq(colour = "#253954", show.legend = FALSE) +
  geom_gene(aes(gene_fill = accver), position = position_strand(.10), width = .045) +
  geom_seq_label(size = 4) +
  coord_collinear() +
  scale_gene_fill_manual(values = first_colours, name = "Sequence") +
  scale_seq_position_continuous(n.breaks = 3, minor_breaks = NULL) +
  theme_ggchord_publication(axis.text = element_text(size = 11)) +
  theme(legend.position = "bottom")
p_first
# ---- firstsave ----
ggsave("first-comparison.png", p_first, width = 10, height = 5, units = "in", dpi = 144)
# ---- end ----
