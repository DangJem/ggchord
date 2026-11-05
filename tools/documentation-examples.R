# Shared executable documentation snippets, canonical website and vignette source.
# Prefer source-backed bundled data; retain metadata outside plotting inputs.
# ---- setup ----
library(dplyr)
seqs <- seq_data_example |>
  select(accver, length)
links <- ribbon_data_example |>
  select(qaccver, saccver, qstart, qend, sstart, send)
# Recover biological strands from the source annotation field.
genes <- gene_data_example |>
  transmute(accver, start, end, strand = source_strand)
# ---- sequencepalette ----
# Optional appearance choices; none of these columns or scales is required input.
seq_palette <- setNames(
  c("#E64B35", "#233E4A", "#269B82", "#D59B17"),
  c("MT108731.1", "MT118296.1", "OQ646790.1", "OR222515.1")
)
gene_palette <- c("+" = "#CE1536", "-" = "#346B9D", "." = "#8291A6")
# ---- featurepalette ----
feature_palette <- c(
  CDS = "#993366", promoter = "#31849B", rep_origin = "#F1E51D",
  primer_bind = "#A020F0", misc_feature = "#A6ACB3", protein_bind = "#31849B"
)
# ---- docstyle ----
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
p_mapped <- ggchord(seq_data_example) +
  geom_seq() +
  geom_link_ribbon(aes(ribbon_fill = pident), data = head(ribbon_data_example, 3)) +
  coord_chord()
p_mapped
# ---- labels ----
label_data <- head(subset(gene_data_example, accver == "OQ646790.1"), 4)
label_data$strand <- label_data$source_strand
p_labels <- ggchord(seq_data_example) +
  geom_seq() +
  geom_gene(data = label_data) +
  geom_gene_label_repel(data = label_data) +
  coord_chord()
p_labels
# ---- trackdata ----
plasmid <- plasmid_example_pUC19
track_seq <- plasmid[c("accver", "length")]
starts <- seq.int(1, plasmid$length, by = 100)
bins <- data.frame(accver = plasmid$accver, start = starts,
                   end = pmin(starts + 99, plasmid$length))
bins$gc <- vapply(seq_len(nrow(bins)), function(i) {
  bases <- strsplit(substr(plasmid$sequence, bins$start[i], bins$end[i]), "")[[1]]
  mean(bases %in% c("G", "C"))
}, numeric(1))
head(bins)
# ---- tracks ----
p_tracks <- ggchord(track_seq) +
  geom_seq() +
  geom_track_bar(aes(track_value = gc), data = bins,
                 position = position_track("gc", width = grid::unit(10, "mm"))) +
  scale_track_value_continuous("gc", limits = c(0, 1), unit = "GC fraction") +
  coord_circular()
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
homolog_source <- bind_rows(plasmid_example_pUC19, plasmid_example_pETDuet_1)
homolog_seq <- select(homolog_source, accver, length)
cluster_features <- homolog_source |>
  find_common_features() |>
  filter(type == "CDS", anno == "AmpR") |>
  mutate(feature_id = paste0(accver, "_AmpR"))
select(cluster_features, feature_id, start, end, strand)
# ---- homologdna ----
reverse_complement <- function(dna) {
  bases <- rev(strsplit(dna, "")[[1]])
  chartr("ACGT", "TGCA", paste(bases, collapse = ""))
}
coding <- left_join(cluster_features, homolog_source, by = "accver") |>
  rowwise() |>
  mutate(dna = substr(sequence, start, end),
         dna = if (strand == "-") reverse_complement(dna) else dna) |>
  ungroup()
cds_dna <- coding$dna
data.frame(feature_id = coding$feature_id, length_bp = nchar(cds_dna),
           identical_to_first = cds_dna == cds_dna[1])
# ---- homologmembership ----
membership <- tibble(
  feature_id = cluster_features$feature_id,
  cluster_id = match(cds_dna, unique(cds_dna))
)
homologs <- link_ggchord_clusters(cluster_features, membership,
  pairs = tibble(qaccver = "pUC19", saccver = "pETDuet-1")
)
membership
# ---- homologs ----
local_homologs <- focus_ggchord_data(homolog_seq, homologs$data, cluster_features,
  loci = select(cluster_features, accver, start, end), expand = 200
)
p_homologs <- ggchord(local_homologs$seq_data) +
  geom_seq() +
  geom_link_ribbon(data = local_homologs$ribbon_data, fill = "grey75") +
  geom_gene(data = local_homologs$gene_data) +
  geom_seq_break() +
  coord_collinear() +
  scale_seq_position_continuous(coordinates = "source")
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
focus_loci <- data.frame(accver = seq_data_example$accver, start = 10000, end = pmin(40000, seq_data_example$length))
focused <- focus_ggchord_data(seq_data_example, ribbon_data_example, loci = focus_loci)
p_focus <- ggchord(focused$seq_data) +
  geom_seq() +
  geom_link_ribbon(data = focused$ribbon_data) +
  geom_seq_break() +
  coord_collinear() +
  scale_seq_position_continuous(coordinates = "source")
p_focus
# ---- segment ----
windows <- data.frame(
  accver = "MT108731.1", start = c(1, 40000), end = c(20000, 55000),
  locus_id = c("left", "right"), orientation = c(1, -1)
)
segmented <- segment_ggchord_data(seq_data_example[1, ], windows = windows)
p_segments <- ggchord(segmented$seq_data) +
  geom_seq() +
  coord_collinear(segments = segmented$segments)
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
seq_real <- read_fasta_lengths("examples/fasta/*.fna")
feature_real <- read_gff3("examples/gff3/*.gff3", unstranded = "keep")
link_real <- read_blast("examples/blastn-standard/*.tsv")
validate_ggchord_data(seq_real, link_real)$summary
# ---- realprepare ----
real_filtered <- filter_ggchord_ribbons(link_real, min_pident = 80, min_length = 500)
real_filtered$report
# ---- realfeatures ----
real_display_features <- feature_real |>
  filter(type == "CDS") |>
  group_by(accver) |>
  slice_max(end - start + 1, n = 6, with_ties = FALSE) |>
  ungroup()
head(select(real_display_features, accver, start, end, strand), 3)
# ---- realchord ----
p_real_chord <- ggchord(seq_real) +
  geom_seq() +
  geom_seq_label() +
  geom_link_ribbon(data = real_filtered$data) +
  geom_feature(data = real_display_features) +
  coord_chord()
p_real_chord
# ---- reallinear ----
p_real_linear <- ggchord(seq_real) +
  geom_seq() +
  geom_seq_label() +
  geom_link_ribbon(data = real_filtered$data) +
  geom_feature(data = feature_real) +
  geom_seq_scale(length = 10000) +
  coord_collinear()
p_real_linear
# ---- plasmidinput ----
library(dplyr)
plasmid_feature_source <- find_common_features(plasmid_example_pUC19)
plasmid_features <- filter(plasmid_feature_source, type != "primer_bind")
plasmid_primers <- find_primer_bindings(plasmid_example_pUC19, set = "universal")
plasmid_sites <- find_restriction_sites(plasmid_example_pUC19)
display_sites <- filter_restriction_sites(plasmid_sites, enzymes = c("EcoRI", "BamHI", "HindIII"))
# ---- plasmid ----
plasmid_position <- position_feature_stack(base_position = position_plasmid())
p_plasmid <- ggchord(plasmid_example_pUC19) +
  geom_seq() +
  geom_seq_center_label() +
  geom_feature_plasmid(data = plasmid_features, position = plasmid_position) +
  geom_feature_label_repel(data = plasmid_features, position = plasmid_position) +
  geom_primer(data = plasmid_primers) +
  geom_restriction_site(data = display_sites) +
  coord_circular()
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
bed_path <- tempfile(fileext = ".bed")
writeLines("pUC19\t866\t1455", bed_path)
bed_features <- read_bed(bed_path)
unlink(bed_path)
bed_features
# ---- formatpaf ----
paf_source <- read_paf("examples/format-inputs/emales-pair.paf")
paf_links <- filter(paf_source, !is.na(mapq), mapq >= 30)
paf_seqs <- bind_rows(
  transmute(paf_source, accver = qaccver, length = qlen),
  transmute(paf_source, accver = saccver, length = slen)
) |> distinct()
head(select(paf_links, qaccver, qstart, qend, saccver, sstart, send), 3)
# ---- drawpaf ----
p_paf <- ggchord(paf_seqs) +
  geom_seq() +
  geom_seq_label() +
  geom_link_ribbon(data = paf_links) +
  coord_collinear()
p_paf
# ---- standardblast ----
blast_source <- read_blast("examples/blastn-standard/MT108731.1__MT118296.1.tsv")
head(blast_source[1:12], 3)
# ---- blastendpoints ----
blast_links <- blast_source[c("qaccver", "saccver", "qstart", "qend", "sstart", "send")]
head(blast_links, 3)
# ---- nativeproteininput ----
protein_seq_data_example[c("accver", "length", "sequence_unit")]
# ---- nativeproteinlinks ----
protein_ribbon_data_example[c("qaccver", "saccver", "qstart", "qend", "sstart", "send", "pident")]
# ---- nativeproteinplot ----
p_native_protein <- ggchord(protein_seq_data_example, protein_ribbon_data_example) +
  geom_seq() +
  geom_seq_label() +
  geom_link_ribbon() +
  geom_seq_scale(length = 50) +
  coord_collinear()
p_native_protein
# ---- nativeproteinchord ----
p_native_chord <- ggchord(protein_seq_data_example, protein_ribbon_data_example) +
  geom_seq() +
  geom_seq_label() +
  geom_link_ribbon() +
  coord_chord()
p_native_chord
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
select(selected_hsp, qaccver, qstart, qend, length, evalue)
# ---- proteinproject ----
query_interval <- tibble(
  feature_id = "query_HSP", parent_id = selected_hsp$qaccver,
  start = selected_hsp$qstart, end = selected_hsp$qend
)
projected_query <- project_ggchord_features(query_interval, emale_cds)
head(projected_query[c("accver", "start", "end", "strand")])
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
# ---- ribboninputs ----
ribbon_sequences <- seq_data_example
observed_links <- ribbon_data_example
# ---- ribbons ----
selected_links <- filter_ggchord_ribbons(observed_links, min_pident = 80, min_length = 500)
unique_links <- deduplicate_ggchord_ribbons(selected_links$data, by = "overlap",
                                           min_reciprocal_overlap = 0.9)
merged_links <- merge_ggchord_ribbons(unique_links$data, max_gap = 100,
                                      min_pident_difference = 2)
bundled_links <- bundle_ggchord_ribbons(observed_links, ribbon_sequences, bins = 30)
# ---- ribbondensity ----
p_density <- ggchord(ribbon_sequences, observed_links) +
  geom_seq() +
  stat_ribbon_density(aes(ribbon_alpha = after_stat(bundle_n)), bins = 4,
                      weight = "count", fill = "#007D74") +
  scale_ribbon_alpha_continuous(name = "Local alignment count", range = c(0.18, 0.95)) +
  coord_chord()
p_bundle <- ggchord(ribbon_sequences, observed_links) +
  geom_seq() +
  stat_ribbon_bundle(bins = 30) +
  coord_chord()
p_density
# ---- pointlinks ----
points <- transform(head(observed_links, 3),
  qpos = floor((qstart + qend) / 2),
  spos = floor((sstart + send) / 2)
)
p_points <- ggchord(seq_data = ribbon_sequences) +
  geom_seq() +
  geom_link_line(data = points, arrow = grid::arrow(length = grid::unit(2, "mm"))) +
  coord_chord()
p_points
# ---- featurestack ----
plasmid <- plasmid_example_pUC19
stack_features <- head(find_common_features(plasmid), 6)
stack_position <- position_feature_stack(base_position = position_plasmid())
p_stack <- ggchord(plasmid) +
  geom_seq() +
  geom_seq_center_label() +
  geom_feature(data = stack_features, position = stack_position) +
  geom_feature_label_repel(data = stack_features, position = stack_position) +
  coord_circular()
p_stack
# ---- labellayouts ----
label_plots <- lapply(c("auto", "radial", "arc"), function(layout) {
  ggchord(seq_data_example) +
    geom_seq() +
    geom_gene(data = label_data) +
    geom_gene_label_repel(data = label_data, label_layout = layout) +
    coord_chord()
})
names(label_plots) <- c("auto", "radial", "arc")
label_plots$auto
# ---- manualgenelabels ----
p_manual_labels <- ggchord(seq_data_example) +
  geom_seq() +
  geom_gene(data = label_data) +
  geom_gene_label(data = label_data, label_rotation = 30) +
  coord_chord()
p_manual_labels
# ---- decorations ----
origin <- subset(find_common_features(plasmid_example_pUC19), type == "rep_origin")
p_decorated <- ggchord(plasmid_example_pUC19) +
  geom_seq() +
  geom_seq_center_label() +
  geom_seq_region(data = origin) +
  coord_circular()
p_decorated
# ---- multitrack ----
line_bins <- transform(bins, pos = floor((start + end) / 2))
line_bins$gc_mean <- as.numeric(stats::filter(bins$gc, rep(1 / 3, 3)))
line_bins <- subset(line_bins, is.finite(gc_mean))
p_multitrack <- p_tracks +
  geom_track_line(aes(track_value = gc_mean), data = line_bins,
                  position = position_track("gc", width = grid::unit(10, "mm")))
focused_bins <- focus_ggchord_track(bins, track_seq,
  loci = data.frame(accver = track_seq$accver, start = 500, end = 1500), type = "interval")
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
source_mark <- map_ggchord_position(p_focus, "MT108731.1", 10100, source = TRUE)
p_marked <- p_focus +
  ggplot2::geom_point(data = source_mark, aes(x, y), inherit.aes = FALSE)
p_marked
# ---- styling ----
p_styled <- ggchord(seq_data_example) +
  geom_seq() +
  geom_link_ribbon(aes(ribbon_fill = qaccver), data = head(ribbon_data_example, 1)) +
  scale_ribbon_fill_manual(values = c("MT108731.1" = "#269B82"), name = "Query") +
  coord_chord()
p_styled
# ---- layoutresults ----
full_layout <- get_chord_layout(p_collinear)
layer_tables <- export_ggchord_layout(p_collinear, include = c("seq", "gene"), original_data = TRUE)
layout_choice <- optimize_ggchord_layout(seqs, links, optimize = c("order", "orientation"), weight = "count")
# ---- parameters ----
p_parameters <- ggchord(seq_data_example) +
  geom_seq(seq_radius = c("MT118296.1" = 0.72)) +
  geom_seq_label() +
  coord_chord()
p_parameters
# ---- parameteruniform ----
p_parameter_uniform <- ggchord(seq_data_example) +
  geom_seq(seq_gap = 0.08) +
  geom_seq_label() +
  coord_chord()
p_parameter_uniform
# ---- parametervector ----
p_parameter_vector <- ggchord(seq_data_example) +
  geom_seq(seq_radius = c(1, 0.8, 1.2, 1)) +
  geom_seq_label() +
  coord_chord()
p_parameter_vector
# ---- parametergenes ----
parameter_genes <- subset(gene_data_example,
  (accver == "OQ646790.1" & start == 43765) |
  (accver == "OR222515.1" & start == 19246))
parameter_genes$strand <- parameter_genes$source_strand
parameter_genes[c("accver", "start", "end", "strand", "anno")]
# ---- parameterstrand ----
p_parameter_strand <- ggchord(seq_data_example) +
  geom_seq() +
  geom_gene(data = parameter_genes,
    position = position_strand(c("+" = 0.08, "-" = -0.08))) +
  coord_chord()
p_parameter_strand
# ---- parameterwidth ----
p_parameter_width <- ggchord(seq_data_example) +
  geom_seq() +
  geom_gene(data = parameter_genes, width = c("+" = 0.03, "-" = 0.08)) +
  coord_chord()
p_parameter_width
# ---- parameterlabels ----
p_parameter_labels <- ggchord(seq_data_example) +
  geom_seq() +
  geom_seq_label(seq_label_radius = c("MT108731.1" = 0.75, "OQ646790.1" = 1.15)) +
  coord_chord()
p_parameter_labels
# ---- frontseqtable ----
seq_data_example
# ---- frontlinks ----
chord_links <- head(ribbon_data_example, 1)
# ---- frontlinktable ----
head(ribbon_data_example[c("qaccver", "qstart", "qend", "saccver", "sstart", "send", "pident")], 5)
# ---- frontgenes ----
chord_genes_example <- subset(gene_data_example, accver == "OQ646790.1" & start == 43765)
chord_genes_example$strand <- chord_genes_example$source_strand
# ---- frontgenetable ----
example_cds <- head(gene_data_example, 5)
example_cds$strand <- example_cds$source_strand
example_cds[c("accver", "start", "end", "strand", "anno")]
# ---- frontselectedgene ----
chord_genes_example[c("accver", "start", "end", "strand", "anno")]
# ---- frontbackbone ----
p_chord_backbone <- ggchord(seq_data = seq_data_example) +
  geom_seq() +
  geom_seq_label() +
  coord_chord()
p_chord_backbone
# ---- frontribbon ----
p_chord_ribbon <- p_chord_backbone +
  geom_link_ribbon(data = chord_links)
p_chord_ribbon
# ---- frontfeature ----
p_chord_feature <- p_chord_ribbon +
  geom_gene(data = chord_genes_example)
p_chord_feature
# ---- frontrotation ----
p_chord_rotated <- p_chord_backbone + coord_chord(rotation = 90)
p_chord_rotated
# ---- frontgap ----
p_chord_gap <- ggchord(seq_data = seq_data_example) +
  geom_seq(seq_gap = c("MT108731.1" = 0.12)) +
  geom_seq_label() +
  coord_chord()
p_chord_gap
# ---- frontradius ----
chord_radius <- c("MT108731.1" = 1, "MT118296.1" = 0.72,
                  "OQ646790.1" = 1, "OR222515.1" = 1)
p_chord_radius <- ggchord(seq_data = seq_data_example) +
  geom_seq(seq_radius = chord_radius) +
  geom_seq_label() +
  coord_chord()
p_chord_radius
# ---- frontcurvature ----
chord_curvature <- c("MT108731.1" = 1, "MT118296.1" = 0,
                     "OQ646790.1" = 1, "OR222515.1" = 1)
p_chord_curvature <- ggchord(seq_data = seq_data_example) +
  geom_seq(seq_curvature = chord_curvature) +
  geom_seq_label() +
  coord_chord()
p_chord_curvature
# ---- frontorientation ----
chord_orientation <- c("MT108731.1" = 1, "MT118296.1" = -1,
                       "OQ646790.1" = 1, "OR222515.1" = 1)
p_chord_orientation <- ggchord(seq_data = seq_data_example) +
  geom_seq(seq_orientation = chord_orientation) +
  geom_seq_label() +
  coord_chord()
p_chord_orientation
# ---- frontroute ----
p_chord_route <- p_chord_backbone +
  geom_link_ribbon(data = chord_links, link_ctrl_point = c(0, 0.8))
p_chord_route
# ---- frontobstacle ----
chord_obstacle <- subset(gene_data_example, accver == "MT108731.1" & start == 48053)
chord_obstacle$strand <- chord_obstacle$source_strand
p_chord_obstacle <- p_chord_backbone +
  geom_gene(data = chord_obstacle, position = position_strand(c("+" = -0.05)))
# ---- frontnoavoid ----
p_chord_no_avoid <- p_chord_obstacle +
  geom_link_ribbon(data = chord_links, link_avoid = "none")
p_chord_no_avoid
# ---- frontsmooth ----
p_chord_smooth <- p_chord_obstacle +
  geom_link_ribbon(data = chord_links, link_avoid = "smooth")
p_chord_smooth
# ---- frontseqlabels ----
p_chord_sequence_labels <- ggchord(seq_data = seq_data_example) +
  geom_seq() +
  geom_seq_label(seq_label_radius = 0.7, seq_label_orientation = "horizontal") +
  coord_chord()
p_chord_sequence_labels
# ---- frontgenelabels ----
p_chord_gene_labels <- ggchord(seq_data_example) +
  geom_seq() +
  geom_seq_label() +
  geom_gene(data = chord_genes_example) +
  geom_gene_label(data = chord_genes_example, label_radial_offset = 0.10) +
  coord_chord()
p_chord_gene_labels
# ---- frontmarkers ----
chord_query_region <- data.frame(
  accver = chord_links$qaccver, start = chord_links$qstart, end = chord_links$qend
)
chord_midpoints <- transform(chord_links,
  qpos = floor((qstart + qend) / 2), spos = floor((sstart + send) / 2)
)
p_chord_markers <- p_chord_backbone +
  geom_seq_region(data = chord_query_region, region_side = "outside") +
  geom_link_line(data = chord_midpoints, link_type = "straight",
                 arrow = grid::arrow(length = grid::unit(2, "mm")))
p_chord_markers
# ---- frontorder ----
chord_order <- c("MT108731.1", "OQ646790.1", "MT118296.1", "OR222515.1")
p_chord_order <- ggchord(seq_data = seq_data_example) +
  geom_seq(seq_order = chord_order) +
  geom_seq_label() +
  coord_chord()
p_chord_order
# ---- frontlinearbase ----
p_linear_base <- ggchord(seq_data = seq_data_example) +
  geom_seq() +
  geom_seq_label() +
  coord_collinear()
p_linear_base
# ---- frontlinearlink ----
p_linear_link <- p_linear_base + geom_link_ribbon(data = chord_links)
p_linear_link
# ---- frontlinearcds ----
p_linear_cds <- p_linear_link + geom_gene(data = chord_genes_example)
p_linear_cds
# ---- frontlineargap ----
p_linear_gap <- ggchord(seq_data = seq_data_example) +
  geom_seq(seq_gap = c("MT108731.1" = 0.45)) +
  geom_seq_label() +
  coord_collinear() +
  geom_link_ribbon(data = chord_links)
p_linear_gap
# ---- frontlinearscale ----
p_linear_scale <- p_linear_link + geom_seq_scale(length = 10000)
p_linear_scale
# ---- frontlinearcross ----
cross_link <- head(subset(ribbon_data_example,
  qaccver == "MT108731.1" & saccver == "OQ646790.1"), 1)
p_linear_cross <- p_linear_link +
  geom_link_ribbon(data = cross_link, adjacent_only = FALSE)
p_linear_cross
# ---- frontcrosswidth ----
p_cross_width <- p_linear_link +
  geom_link_ribbon(data = cross_link, adjacent_only = FALSE,
    nonadjacent_linewidth = 0.8)
p_cross_width
# ---- frontcrosssame ----
p_cross_same <- p_linear_link +
  geom_link_ribbon(data = cross_link, adjacent_only = FALSE,
    linewidth = 0.2, nonadjacent_linewidth = 0)
p_cross_same
# ---- frontimportseq ----
head(read_fasta_lengths(files = "examples/fasta/*.fna"), 4)
# ---- frontimportgff ----
gff <- read_gff3(files = "examples/gff3/*.gff3", unstranded = "keep")
head(gff[c("accver", "start", "end", "strand", "type")], 4)
# ---- frontimportblast ----
hits <- read_blast(file = "examples/blastn-standard/MT108731.1__MT118296.1.tsv")
head(hits[c("qaccver", "saccver", "qstart", "qend", "sstart", "send")], 3)
# ---- frontimportgr ----
gr <- GenomicRanges::GRanges("a", IRanges::IRanges(10, 30), strand = "+", seqlengths = c(a = 100))
as_ggchord_feature(gr)
# ---- frontimportgi ----
query <- GenomicRanges::GRanges("a", IRanges::IRanges(10, 30), seqlengths = c(a = 100))
subject <- GenomicRanges::GRanges("b", IRanges::IRanges(40, 60), seqlengths = c(b = 120))
pair <- InteractionSet::GInteractions(query, subject)
as_ggchord_link(pair)
# ---- frontpreview ----
view_ggchord(p_linear_base, width = 8, height = 5, units = "in")
# ---- frontsave ----
ggsave("comparison.png", p_linear_base, width = 8, height = 5, units = "in", dpi = 144)
# ---- frontgeometry ----
layout <- export_ggchord_layout(p_linear_base, include = "seq")
head(unique(layout$seq[c("layer_id", "component", "source_row")]), 4)
# ---- firstcomparison ----
lesson_ids <- c("MT108731.1", "MT118296.1")
lesson_sequences <- subset(seq_data_example, accver %in% lesson_ids)
lesson_links <- subset(ribbon_data_example, qaccver %in% lesson_ids & saccver %in% lesson_ids)
# ---- firstchord ----
p_first_chord <- ggchord(seq_data = lesson_sequences, ribbon_data = lesson_links) +
  geom_seq() +
  geom_link_ribbon() +
  geom_seq_label() +
  coord_chord()
p_first_chord
# ---- firstlinear ----
p_first_linear <- ggchord(seq_data = lesson_sequences, ribbon_data = lesson_links) +
  geom_seq() +
  geom_link_ribbon() +
  geom_seq_label() +
  coord_collinear()
p_first_linear
# ---- firstcircle ----
plasmid <- plasmid_example_pUC19
p_first_circle <- ggchord(seq_data = plasmid) +
  geom_seq() +
  geom_seq_center_label() +
  coord_circular()
p_first_circle
# ---- firstfeatures ----
features <- find_common_features(plasmid)
origins <- subset(features, type == "rep_origin")
p_first_features <- p_first_circle +
  geom_feature_plasmid(data = origins)
p_first_features
# ---- frontcirclestyle ----
p_circle_single <- ggchord(seq_data = plasmid) +
  geom_seq(seq_style = "single") +
  geom_seq_center_label() +
  coord_circular()
p_circle_single
# ---- frontcirclerotation ----
p_circle_rotated <- p_first_features + coord_circular(rotation = 180)
p_circle_rotated
# ---- practicegenes ----
# Restore the biological strand saved separately in this bundled teaching table.
lesson_genes <- subset(gene_data_example, accver %in% lesson_ids)
lesson_genes$strand <- lesson_genes$source_strand
# Select one CDS per sequence. Prefer a specific annotation when available.
lesson_genes <- subset(lesson_genes, !is.na(anno))
lesson_genes <- lesson_genes[order(
  lesson_genes$anno == "hypothetical protein",
  -(lesson_genes$end - lesson_genes$start)
), ]
lesson_genes <- lesson_genes[!duplicated(lesson_genes$accver), ]
lesson_position <- position_strand(0.09)
# ---- practicechord ----
p_annotated_chord <- p_first_chord +
  geom_gene(data = lesson_genes, position = lesson_position) +
  geom_gene_label(data = lesson_genes, position = lesson_position)
p_annotated_chord
# ---- practicelinear ----
p_annotated_linear <- p_first_linear +
  geom_gene(data = lesson_genes, position = lesson_position) +
  geom_gene_label(data = lesson_genes, position = lesson_position) +
  geom_seq_scale(length = 10000)
p_annotated_linear
# ---- inspectlinks ----
head(lesson_links[c("qaccver", "qstart", "qend", "saccver", "sstart", "send")], 2)
# ---- inspectgenes ----
lesson_genes[c("accver", "start", "end", "strand", "anno")]
# ---- inspectorigin ----
origins[c("accver", "start", "end", "type", "anno")]
# ---- practiceplasmid ----
p_annotated_plasmid <- p_first_features +
  geom_feature_label(
    data = origins, label_orientation = "horizontal",
    label_side = "inside", label_radial_offset = 0.12, size = 4
  )
p_annotated_plasmid
# ---- firstsavebasic ----
ggsave("first-comparison.png", p_first_linear, width = 8, height = 5, units = "in")
# ---- firstdata ----
library(dplyr)
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
  dplyr::slice(c(1, n())) |>
  ungroup()
show_features <- bind_rows(show_features, show_reverse) |>
  distinct(accver, start, end, strand, .keep_all = TRUE) |>
  arrange(accver, start)
show_labels <- show_features
show_optimization <- optimize_ggchord_layout(
  show_sequences, show_links,
  optimize = "orientation", weight = "length"
)
# ---- chordshowcasegenes ----
# Sample each observed strand separately across the source coordinates.
# Keep short CDS and retain the original strand and coordinates.
# Prefer a CDS overlapping a displayed ribbon endpoint within each positive-strand bin.
chord_link_ends <- bind_rows(
  transmute(show_links, accver = qaccver, lo = pmin(qstart, qend), hi = pmax(qstart, qend)),
  transmute(show_links, accver = saccver, lo = pmin(sstart, send), hi = pmax(sstart, send))
)
chord_endpoint_genes <- show_genes |>
  filter(strand == "+") |>
  semi_join(chord_link_ends, by = join_by(accver, overlaps(start, end, lo, hi)))
chord_genes <- show_genes |>
  filter(strand %in% c("+", "-")) |>
  left_join(mutate(chord_endpoint_genes, at_endpoint = TRUE),
    by = c("accver", "start", "end", "strand", "anno")) |>
  group_by(accver, strand) |>
  arrange(start, .by_group = TRUE) |>
  mutate(display_bin = ntile(row_number(), if (first(strand) == "+") 4L else 3L)) |>
  group_by(accver, strand, display_bin) |>
  arrange(desc(coalesce(at_endpoint, FALSE)), is.na(anno) | anno == "hypothetical protein",
    start, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  select(accver, start, end, strand, anno)
chord_labels <- chord_genes
chord_position <- position_strand(c("+" = -0.030, "-" = 0.030))
# ---- showcasechord ----
show_curvature <- c(
  "MT108731.1" = 0, "MT118296.1" = 1,
  "OQ646790.1" = -0.65, "OR222515.1" = 0.75
)
show_radius <- c(
  "MT108731.1" = 1.08, "MT118296.1" = 0.82,
  "OQ646790.1" = 1.18, "OR222515.1" = 1.02
)
p_show_chord <- ggchord(seq_data = show_sequences, ribbon_data = show_links) +
  geom_link_ribbon(link_avoid = "smooth", alpha = 0.46, linewidth = 0.14) +
  geom_seq(
    seq_orientation = show_optimization$seq_orientation,
    seq_curvature = show_curvature, seq_radius = show_radius,
    seq_gap = 0.06, linewidth = 1.1
  ) +
  geom_gene(data = chord_genes, position = chord_position, width = 0.025) +
  geom_gene_label_repel(
    data = chord_labels, label_layout = "radial", label_wrap = 24,
    size = 2.8, label_max_lines = 3, position = chord_position
  ) +
  geom_seq_label(size = 4, seq_label_radius = 0.77, seq_label_hjust = 0.5) +
  scale_seq_colour_manual(
    values = seq_palette, name = "Sequence",
    guide = sequence_guide
  ) +
  scale_gene_fill_manual(values = gene_palette, breaks = c("+", "-"), name = "Strand") +
  guides(gene_fill = strand_guide, ribbon_fill = identity_guide) +
  scale_seq_position_continuous(breaks = c(0, 20000, 40000), minor_breaks = NULL) +
  coord_chord(rotation = 50) +
  doc_theme +
  labs(
    title = "Four genomes, one connected view",
    subtitle = "Mixed backbone geometry · CDS sampled by source strand"
  )
p_show_chord
# ---- showcaseplasmiddata ----
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
  coord_circular(rotation = 90) +
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
  guides(gene_fill = strand_guide, ribbon_fill = identity_guide) +
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
  sites <- find_restriction_sites(sequence)
  # Pick the nearest observed cut to each of three evenly spaced targets.
  targets <- sequence$length * c(1, 3, 5) / 6
  nearest <- vapply(targets, function(x) which.min(abs(sites$position - x)), integer(1))
  sites <- sites[unique(nearest), ]
  ggchord(seq_data = sequence |> select(accver, length)) +
    geom_seq(seq_style = "double", colour = "#252525") +
    geom_feature_plasmid(data = features, position = placement) +
    geom_feature_label_repel(
      data = features, position = placement, external = TRUE,
      label_fit = "shrink", label_min_size = 2.5
    ) +
    geom_restriction_site(data = sites, label_side = "outside", label_size = 2.5) +
    geom_seq_center_label() +
    coord_circular(rotation = 90) +
    theme_ggchord_plasmid() +
    theme(
      panel.background = element_blank(),
      plot.background = element_rect(fill = doc_background, colour = NA)
    )
}
# ---- catalogbatch ----
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
# ---- frontplasmidsummary ----
data.frame(features = nrow(plasmid_features), primers = nrow(plasmid_primers),
           displayed_sites = nrow(display_sites))
# ---- frontduetsummary ----
data.frame(features = nrow(show_plasmid_features), primers = nrow(show_plasmid_primers),
           selected_site_rows = nrow(show_plasmid_sites))
# ---- frontgctable ----
head(bins, 4)
# ---- frontlayoutchoice ----
data.frame(accver = names(layout_choice$seq_orientation),
           orientation = unname(layout_choice$seq_orientation))
# ---- end ----
