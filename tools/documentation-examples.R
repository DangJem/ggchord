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
  c("#E64B35", "#233E4A", "#269B82", "#D59B17"),
  c("MT108731.1", "MT118296.1", "OQ646790.1", "OR222515.1")
)
gene_palette <- c("+" = "#CE1536", "-" = "#346B9D", "." = "#8291A6")
# ---- featurepalette ----
library(ggchord)
feature_palette <- c(
  CDS = "#993366", promoter = "#31849B", rep_origin = "#F1E51D",
  primer_bind = "#A020F0", misc_feature = "#A6ACB3", protein_bind = "#31849B"
)
# ---- docstyle ----
library(ggchord)
doc_background <- "#FFFCF5"
identity_palette <- c("#303B71", "#36537D", "#269B82", "#A5D52C", "#ECEB31")
identity_values <- c(0, 0.5, 0.8, 0.9, 1)
identity_breaks <- c(0, 50, 80, 90, 100)
doc_legend_theme <- theme(
  legend.title = element_text(size = 13, face = "bold"),
  legend.text = element_text(size = 11),
  legend.key.width = unit(9, "mm"), legend.key.height = unit(5, "mm"),
  legend.spacing.x = unit(9, "mm"), legend.spacing.y = unit(7, "mm"),
  legend.margin = margin(l = 4, r = 4, unit = "mm")
)
sequence_guide <- guide_ggchord_legend(nrow = 1, theme = doc_legend_theme)
strand_guide <- guide_ggchord_legend(
  nrow = 1, override.aes = list(strand = c("+", "-")), theme = doc_legend_theme
)
identity_guide <- guide_ggchord_colourbar(
  position = "left", direction = "vertical",
  theme = theme(
    legend.key.height = unit(75, "mm"), legend.key.width = unit(5, "mm"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11), legend.margin = margin(r = 7, unit = "mm")
  )
)
doc_theme <- theme_ggchord_publication(
  axis.text = element_text(size = 11),
  axis.text.check.overlap = TRUE
) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.width = unit(9, "mm"), legend.key.height = unit(5, "mm"),
    legend.spacing.x = unit(9, "mm"), legend.spacing.y = unit(7, "mm"),
    legend.box.spacing = unit(6, "mm"),
    plot.background = element_rect(fill = doc_background, colour = NA),
    panel.background = element_blank(),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 10, colour = "#5C625A"),
    plot.margin = margin(14, 18, 14, 18)
  )
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
  doc_theme
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
  scale_ribbon_fill_gradientn(
    colours = identity_palette, values = identity_values,
    limits = c(0, 100), breaks = identity_breaks, name = "Identity (%)", guide = identity_guide
  ) +
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
  scale_track_fill_gradientn("gc", colours = c("#E4F3EA", "#65B79B", "#27775C")) +
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
  geom_primer(data = primers, colour = "#A020F0") +
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
  geom_link_ribbon(fill = "#269B82", alpha = .65) +
  geom_seq(colour = "#233E4A", show.legend = FALSE) +
  geom_feature(
    data = local_homologs$gene_data, feature_fill = "#F1D334",
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
  geom_link_ribbon(fill = "#6BAFA6", alpha = .5) +
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
  scale_ribbon_fill_gradientn(
    colours = identity_palette, values = identity_values,
    limits = c(0, 100), breaks = identity_breaks, name = "Identity (%)", guide = identity_guide
  ) +
  doc_theme
p_real_chord
# ---- reallinear ----
p_real_linear <- ggchord(seq_real, real_filtered$data, validate = "none") +
  geom_link_ribbon(fill = "#6BAFA6", alpha = .45) +
  geom_seq() +
  geom_seq_label(size = 4) +
  geom_feature(data = feature_real, feature_fill = "#CE1536", show.legend = FALSE, position = position_strand(.18)) +
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
  geom_primer(data = plasmid_primers, colour = "#A020F0") +
  geom_primer_label_repel(data = plasmid_primers, colour = "#A020F0", max_overlaps = 0, size = 3) +
  geom_restriction_site(data = display_sites, colour = "#29364B", label_size = 3, leader = "radial") +
  geom_seq_center_label() +
  scale_feature_fill_manual(values = feature_palette) +
  scale_seq_position_continuous() +
  coord_circular(rotation = 90) +
  theme_ggchord_plasmid() +
  theme(
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.spacing.x = unit(9, "mm"), legend.spacing.y = unit(7, "mm"),
    plot.background = element_rect(fill = doc_background, colour = NA),
    panel.background = element_blank()
  )
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
  geom_link_ribbon(fill = "#6BAFA6", alpha = .5) +
  geom_seq(colour = "#233E4A") +
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
  stat_ribbon_density(aes(ribbon_alpha = after_stat(density)), fill = "#269B82") +
  scale_ribbon_alpha_continuous(
    name = "Relative density", range = c(.1, .6),
    guide = guide_ggchord_legend(position = "bottom", nrow = 1, theme = doc_legend_theme)
  ) +
  scale_seq_colour_manual(values = seq_palette) +
  coord_chord() +
  doc_theme
p_bundle <- ggchord(seq_real, link_real, validate = "none") +
  geom_seq() +
  stat_ribbon_bundle(bins = 30, fill = "#269B82") +
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
    colours = identity_palette, values = identity_values,
    limits = c(0, 100), breaks = identity_breaks,
    name = "Identity (%)", guide = guide_ggchord_colourbar(
      available_aes = "link_colour",
      position = "bottom", direction = "horizontal",
      theme = doc_legend_theme + theme(legend.key.width = unit(65, "mm"), legend.key.height = unit(5, "mm"))
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
    scale_seq_colour_manual(values = seq_palette, name = "Sequence",
      guide = guide_ggchord_legend(nrow = 2, theme = doc_legend_theme)) +
    scale_gene_fill_manual(values = gene_palette, breaks = c("+", "-"), name = "Strand", guide = strand_guide) +
    coord_chord(xlim = c(-1.8, 1.8), ylim = c(-1.6, 1.6)) +
    theme_ggchord_publication(axis = element_blank()) +
    doc_legend_theme +
    theme(
      legend.position = "bottom",
      plot.background = element_rect(fill = doc_background, colour = NA),
      panel.background = element_blank()
    ) +
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
  geom_link_ribbon(fill = "#6BAFA6", alpha = .4) +
  geom_seq() +
  geom_seq_label(size = 4) +
  scale_region_fill_manual(values = c("alignment interval" = "#F4D35E")) +
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
    position = position_track("gc", width = unit(10, "mm")), colour = "#233E4A"
  ) +
  scale_track_value_continuous("gc", limits = c(0, 1), unit = "GC fraction") +
  scale_track_fill_gradientn("gc", colours = c("#E4F3EA", "#65B79B", "#27775C")) +
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
    inherit.aes = FALSE, colour = "#D59B17", size = 3
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
  guides(gene_fill = guide_ggchord_legend(nrow = 1, theme = doc_legend_theme)) +
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
first_colours <- setNames(c("#E64B35", "#233E4A"), ids)
p_first <- ggchord(
  seq_data = first_sequences, ribbon_data = first_links,
  gene_data = first_genes
) +
  geom_link_ribbon(fill = "#82B6AD", alpha = .55) +
  geom_seq(colour = "#233E4A", show.legend = FALSE) +
  geom_gene(aes(gene_fill = accver), position = position_strand(.10), width = .045) +
  geom_seq_label(size = 4) +
  coord_collinear() +
  scale_gene_fill_manual(
    values = first_colours, name = "Sequence",
    guide = guide_ggchord_legend(nrow = 1, theme = theme(
      legend.title = element_text(size = 13, face = "bold"), legend.text = element_text(size = 11),
      legend.key.width = unit(9, "mm"), legend.key.height = unit(5, "mm"), legend.margin = margin(l = 7, unit = "mm")
    ))
  ) +
  scale_seq_position_continuous(n.breaks = 3, minor_breaks = NULL) +
  theme_ggchord_publication(axis.text = element_text(size = 11)) +
  theme(
    legend.position = "bottom",
    plot.background = element_rect(fill = "#FFFCF5", colour = NA),
    panel.background = element_blank()
  )
p_first
# ---- firstsave ----
ggsave("first-comparison.png", p_first, width = 10, height = 5, units = "in", dpi = 144)
# ---- showcasedata ----
library(ggchord)
library(dplyr)
library(purrr)
library(tidyr)
# Original repository records; extract the tutorial input bundle first.
show_sequences <- read_fasta_lengths(files = "examples/fasta/*.fna") |>
  select(accver, length)
show_annotations <- read_gff3(files = "examples/gff3/*.gff3", unstranded = "keep")
show_genes <- show_annotations |>
  filter(type == "CDS") |>
  select(accver, start, end, strand, anno)
show_alignments <- read_blast(files = "examples/blastn-standard/*.tsv")
# At most three substantial observed HSPs per pair; no coordinates are invented.
show_links <- show_alignments |>
  filter(length >= 900, pident >= 75) |>
  group_by(qaccver, saccver) |>
  slice_max(length, n = 3, with_ties = FALSE) |>
  ungroup()
# Select up to four source CDS per sequence, spread across annotation order.
# No feature is moved or invented; keep complete source tables separately.
show_features <- show_genes |>
  filter(end - start + 1 >= 500) |>
  group_by(accver) |>
  arrange(start, .by_group = TRUE) |>
  mutate(display_bin = ntile(row_number(), 4), informative = !is.na(anno) & anno != "hypothetical protein") |>
  group_by(accver, display_bin) |>
  arrange(desc(informative), desc(end - start), .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  select(accver, start, end, strand, anno)
# Also retain well-separated negative-strand CDS, if not already represented.
show_reverse <- show_genes |>
  filter(strand == "-", end - start + 1 >= 500) |>
  group_by(accver) |>
  arrange(start, .by_group = TRUE) |>
  slice(c(1, n())) |>
  ungroup()
show_features <- bind_rows(show_features, show_reverse) |>
  distinct(accver, start, end, strand, .keep_all = TRUE) |>
  arrange(accver, start)
show_labels <- show_features
show_optimization <- optimize_ggchord_layout(
  show_sequences, show_links,
  optimize = "orientation", weight = "length"
)
# ---- showcasechord ----
show_curvature <- c(
  "MT108731.1" = 0, "MT118296.1" = 0.85,
  "OQ646790.1" = 1.1, "OR222515.1" = -0.4
)
show_radius <- c(
  "MT108731.1" = 1, "MT118296.1" = 1.05,
  "OQ646790.1" = 0.95, "OR222515.1" = 1.08
)
p_show_chord <- ggchord(show_sequences, show_links, show_features) +
  geom_link_ribbon(aes(ribbon_fill = pident), link_avoid = "smooth", alpha = 0.58, linewidth = 0.18) +
  geom_seq(
    seq_orientation = show_optimization$seq_orientation,
    seq_curvature = show_curvature, seq_radius = show_radius,
    seq_gap = 0.06, linewidth = 1.1
  ) +
  geom_gene(position = position_strand(0.07), width = 0.035) +
  geom_gene_label_repel(
    data = show_labels, label_layout = "radial", label_wrap = 24,
    size = 3.1, label_max_lines = 3, position = position_strand(0.07)
  ) +
  geom_seq_label(size = 4, seq_label_radius = 0.84, seq_label_hjust = 0.5) +
  scale_seq_colour_manual(
    values = seq_palette, name = "Sequence",
    guide = sequence_guide
  ) +
  scale_gene_fill_manual(values = gene_palette, breaks = c("+", "-"), name = "Strand") +
  scale_ribbon_fill_gradientn(
    colours = identity_palette, values = identity_values,
    limits = c(0, 100), breaks = identity_breaks, name = "Identity (%)", guide = identity_guide
  ) +
  guides(gene_fill = strand_guide) +
  scale_seq_position_continuous(breaks = c(0, 20000, 40000), minor_breaks = NULL) +
  coord_chord(rotation = 50) +
  doc_theme +
  labs(
    title = "Four genomes, one connected view",
    subtitle = "Signed curvature and radius · selected strand-aware CDS · local ribbon clearance"
  )
p_show_chord
# ---- showcaseplasmiddata ----
library(ggchord)
library(dplyr)
data(plasmid_example_pETDuet_1)
show_plasmid <- plasmid_example_pETDuet_1 |>
  select(accver, length)
show_plasmid_features <- find_common_features(plasmid_example_pETDuet_1)
show_plasmid_primers <- find_primer_bindings(plasmid_example_pETDuet_1, set = "reference")
show_plasmid_all_sites <- find_restriction_sites(plasmid_example_pETDuet_1)
show_plasmid_sites <- filter_restriction_sites(
  show_plasmid_all_sites,
  set = "unique_6plus", parent_set = "commercial_all"
)
# ---- showcaseplasmid ----
show_plasmid_position <- position_feature_stack(
  spacing = 0.085, base_position = position_plasmid()
)
p_show_plasmid <- ggchord(show_plasmid, validate = "none") +
  geom_seq(seq_style = "double", colour = "#252525") +
  geom_feature_plasmid(data = show_plasmid_features, position = show_plasmid_position) +
  geom_feature_label_repel(
    data = show_plasmid_features, position = show_plasmid_position,
    external = TRUE, max_overlaps = 0,
    label_fit = "shrink", label_min_size = 2.5
  ) +
  geom_primer(data = show_plasmid_primers, colour = "#A020F0") +
  geom_primer_label_repel(
    data = show_plasmid_primers, colour = "#A020F0", max_overlaps = 0
  ) +
  geom_restriction_site(data = show_plasmid_sites, colour = "#252525") +
  geom_seq_center_label() +
  coord_circular(
    rotation = 90, fit = "manual",
    xlim = c(-3.4, 4.6), ylim = c(-1.6, 2.4)
  ) +
  theme_ggchord_plasmid() +
  theme(
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.spacing.x = unit(9, "mm"), legend.spacing.y = unit(7, "mm"),
    plot.background = element_rect(fill = doc_background, colour = NA),
    panel.background = element_blank()
  )
p_show_plasmid
# ---- showcasegc ----
# Quantitative GC tracks remain in the dedicated tracks tutorial.
# The flagship comparison reserves the space for genes and routed alignments.
# ---- showcasecollinear ----
# Omit two short, crowded candidates from this linear display only.
# Both source rows remain in show_features and in the chord display.
show_linear_features <- show_features |>
  filter(!anno %in% c("tail assembly protein", "gp4 protein"))
# Separate adjacent and cross-row links, retaining original endpoint values.
show_pair_rows <- show_links |>
  mutate(
    query_row = match(qaccver, show_sequences$accver),
    subject_row = match(saccver, show_sequences$accver),
    span = abs(query_row - subject_row)
  )
show_adjacent <- show_pair_rows |> filter(span == 1)
# One longest observed HSP per non-neighbour pair keeps the long routes legible.
show_cross <- show_pair_rows |>
  filter(span > 1) |>
  group_by(qaccver, saccver) |>
  slice_max(length, n = 1, with_ties = FALSE) |>
  ungroup()
# Inspect a backbone-only layout using the same alignment and gap as the final view.
show_linear_coord <- coord_collinear(align = "center", orientation = "auto", ratio = 0.7)
show_routing_base <- ggchord(show_sequences, show_links) +
  geom_seq(seq_gap = 0.42) + show_linear_coord
show_reference <- get_chord_layout(show_routing_base)$sequence_reference
# Shared controls sit just beyond the intervening backbones, rather than far away.
show_cross_controls <- map2(show_cross$query_row, show_cross$subject_row, function(query_row, subject_row) {
  between <- seq(min(query_row, subject_row) + 1L, max(query_row, subject_row) - 1L)
  middle_ids <- show_sequences$accver[between]
  right_edge <- max(show_reference$x_offset[middle_ids] + show_reference$lens[middle_ids] * show_reference$x_scale)
  middle_y <- mean(show_reference$y[show_sequences$accver[c(query_row, subject_row)]])
  c(right_edge + if (abs(query_row - subject_row) == 3) 0.85 else 0.65, middle_y)
})
p_show_collinear <- ggchord(show_sequences, gene_data = show_linear_features) +
  geom_link_ribbon(
    aes(ribbon_fill = pident),
    data = show_adjacent,
    link_avoid = "smooth", alpha = 0.48, linewidth = 0.22
  ) +
  geom_link_ribbon(
    aes(ribbon_fill = pident),
    data = show_cross,
    adjacent_only = FALSE, link_ctrl_point = show_cross_controls,
    link_avoid = "smooth", alpha = 0.28, linewidth = 0.3
  ) +
  geom_seq(seq_gap = 0.42, linewidth = 1) +
  geom_gene(position = position_strand(0.09), width = 0.045) +
  geom_gene_label_repel(
    data = show_linear_features |> arrange(desc(nchar(anno)), accver, start),
    label_layout = "auto", label_wrap = 24,
    label_side = "auto", position = position_strand(0.09), size = 3.3
  ) +
  geom_seq_label(size = 4) +
  geom_seq_scale(length = 10000, offset = unit(c(0, 10), "mm")) +
  scale_seq_colour_manual(
    values = seq_palette, name = "Sequence",
    guide = sequence_guide
  ) +
  scale_gene_fill_manual(values = gene_palette, breaks = c("+", "-"), name = "Strand") +
  scale_ribbon_fill_gradientn(
    colours = identity_palette, values = identity_values,
    limits = c(0, 100), breaks = identity_breaks, name = "Identity (%)", guide = identity_guide
  ) +
  guides(gene_fill = strand_guide) +
  scale_seq_position_continuous(breaks = c(0, 20000, 40000, 60000, 80000), minor_breaks = NULL) +
  show_linear_coord +
  doc_theme +
  labs(
    title = "Sequence architecture, across four genomes",
    subtitle = "Selected labelled CDS · strand-separated arrows · routed cross-row comparisons"
  )
p_show_collinear
# ---- learnseq ----
p_learn_seq <- ggchord(seq_data = show_sequences) +
  geom_seq() +
  geom_seq_label(size = 4) +
  scale_seq_colour_manual(
    values = seq_palette, name = "Sequence",
    guide = sequence_guide
  ) +
  scale_seq_position_continuous(breaks = c(0, 20000, 40000), minor_breaks = NULL) +
  coord_chord() +
  doc_theme
p_learn_seq
# ---- learnlinks ----
p_learn_links <- p_learn_seq +
  geom_link_ribbon(
    aes(ribbon_fill = pident),
    data = show_links,
    adjacent_only = FALSE, link_avoid = "smooth"
  ) +
  scale_ribbon_fill_gradientn(
    colours = identity_palette, values = identity_values, limits = c(0, 100),
    breaks = identity_breaks, name = "Identity (%)", guide = identity_guide
  )
p_learn_links
# ---- learngenes ----
p_learn_genes <- p_learn_links +
  geom_gene(data = show_features, position = position_strand(0.08)) +
  scale_gene_fill_manual(values = gene_palette, breaks = c("+", "-"), name = "Strand") +
  guides(gene_fill = strand_guide)
p_learn_genes
# ---- learnlabels ----
p_learn_labels <- p_learn_genes +
  geom_gene_label_repel(
    data = show_features, label_layout = "arc", label_wrap = 24,
    label_max_lines = 3, position = position_strand(0.08), size = 3
  )
p_learn_labels
# ---- learnlinear ----
# Reuse exactly the same data and layers; replace only the coordinate system.
p_linear_seq <- p_learn_seq + coord_collinear()
p_linear_links <- p_learn_links + coord_collinear()
p_linear_labels <- p_learn_labels + coord_collinear(ratio = 0.65)
p_linear_labels
# ---- learnplasmid ----
plasmid_position <- position_feature_stack(spacing = 0.085, base_position = position_plasmid())
p_plasmid_backbone <- ggchord(seq_data = plasmid_seq) +
  geom_seq(seq_style = "double", colour = "#252525") +
  geom_seq_center_label() +
  coord_circular(rotation = 90) +
  theme_ggchord_plasmid() +
  theme(
    panel.background = element_blank(),
    plot.background = element_rect(fill = doc_background, colour = NA)
  )
p_plasmid_features <- p_plasmid_backbone +
  geom_feature_plasmid(data = plasmid_features, position = plasmid_position)
p_plasmid_text <- p_plasmid_features +
  geom_feature_label_repel(
    data = plasmid_features, position = plasmid_position,
    external = TRUE, label_fit = "shrink", label_min_size = 2.5
  )
p_plasmid_text
# ---- catalogdata ----
library(ggchord)
library(dplyr)
library(purrr)
# Thirteen unique records; the compatibility spelling pUC19c is not an extra plasmid.
catalog_objects <- c(
  "plasmid_example_pSB1C3", "plasmid_example_pUC19", "plasmid_example_pBR322",
  "plasmid_example_pBluescript_II_SK_plus", "plasmid_example_pET_28a_plus",
  "plasmid_example_pETDuet_1", "plasmid_example_pDONR221",
  "plasmid_example_pcDNA3_1_plus", "plasmid_example_pTRE_Tight_BI",
  "plasmid_example_pTRIPZ", "plasmid_example_pCAMBIA1300",
  "plasmid_example_pEarleyGate_201", "plasmid_example_pSpCas9_BB_2A_GFP_PX458"
)
read_catalog_record <- function(object) {
  target <- new.env(parent = emptyenv())
  data(list = object, package = "ggchord", envir = target)
  target[[object]]
}
catalog_sequences <- map(catalog_objects, read_catalog_record)
catalog_features <- map(catalog_sequences, find_common_features)
catalog_summary <- map2(catalog_sequences, catalog_features, function(sequence, features) {
  tibble(accver = sequence$accver, length = sequence$length, features = nrow(features))
}) |>
  bind_rows()
# ---- catalogplots ----
make_catalog_plot <- function(sequence, features) {
  placement <- position_feature_stack(spacing = 0.085, base_position = position_plasmid())
  ggchord(seq_data = sequence |> select(accver, length)) +
    geom_seq(seq_style = "double", colour = "#252525") +
    geom_feature_plasmid(data = features, position = placement) +
    geom_feature_label_repel(
      data = features, position = placement, external = FALSE,
      label_fit = "shrink", label_min_size = 2.5, max_overlaps = 0
    ) +
    geom_seq_center_label() +
    coord_circular(rotation = 90, fit = "manual", xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35)) +
    theme_ggchord_plasmid() +
    theme(
      panel.background = element_blank(),
      plot.background = element_rect(fill = doc_background, colour = NA)
    )
}
catalog_plots <- map2(catalog_sequences, catalog_features, make_catalog_plot)
names(catalog_plots) <- catalog_summary$accver
# Save any record at the same size used in the gallery.
# ggsave("pBR322.png", catalog_plots[["pBR322"]], width = 10, height = 8, dpi = 144)
# ---- optionalinteractions ----
interaction_inputs <- NULL
if (requireNamespace("GenomicRanges", quietly = TRUE) &&
  requireNamespace("InteractionSet", quietly = TRUE)) {
  source_links <- ribbon_data_example |>
    slice_head(n = 2)
  lengths <- setNames(seq_data_example$length, seq_data_example$accver)
  query <- GenomicRanges::GRanges(
    seqnames = source_links$qaccver,
    ranges = IRanges::IRanges(
      pmin(source_links$qstart, source_links$qend),
      pmax(source_links$qstart, source_links$qend)
    ),
    strand = ifelse(source_links$qstart <= source_links$qend, "+", "-"),
    seqlengths = lengths
  )
  subject <- GenomicRanges::GRanges(
    seqnames = source_links$saccver,
    ranges = IRanges::IRanges(
      pmin(source_links$sstart, source_links$send),
      pmax(source_links$sstart, source_links$send)
    ),
    strand = ifelse(source_links$sstart <= source_links$send, "+", "-"),
    seqlengths = lengths
  )
  interactions <- InteractionSet::GInteractions(query, subject)
  interaction_inputs <- list(
    seq_data = as_ggchord_seq(interactions),
    ribbon_data = as_ggchord_link(interactions)
  )
}
# ---- inspectends ----
link_endpoints <- show_links |>
  select(qaccver, saccver, qstart, qend, sstart, send) |>
  pivot_longer(
    cols = c(qstart, qend, sstart, send),
    names_to = c("side", "endpoint"), names_pattern = "([qs])(start|end)",
    values_to = "position"
  )
head(link_endpoints)
# ---- end ----
