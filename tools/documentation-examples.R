# Shared executable documentation snippets, v0.16.5.
# Prefer source-backed bundled data; retain metadata outside plotting inputs.
# ---- setup ----
library(ggchord)
data(seq_data_example)
data(ribbon_data_example)
data(gene_data_example)
seqs <- seq_data_example[c("accver", "length")]
links <- ribbon_data_example[c("qaccver", "saccver", "qstart", "qend", "sstart", "send")]
# The fixture preserves biological strands separately from its old display strands.
genes <- transform(gene_data_example, strand = source_strand)[c("accver", "start", "end", "strand")]
# Optional appearance choices; none of these columns or scales is required input.
seq_palette <- setNames(c("#487A70", "#B48A61", "#637E9D", "#A16E7D"), seqs$accver)
gene_palette <- c("+" = "#56887C", "-" = "#B58169", "." = "#AAB9AE")
feature_palette <- c(CDS = "#6E9A8A", promoter = "#8798B2", rep_origin = "#C4A46A",
  primer_bind = "#A28EAA", misc_feature = "#CDD5CE", protein_bind = "#A9B9C4")
doc_theme <- theme_ggchord_publication(axis.text = element_text(size = 9),
  axis.text.check.overlap = TRUE) + theme(legend.position = "bottom")
doc_axis <- function(coordinates = "relative") scale_seq_position_continuous(
  coordinates = coordinates, n.breaks = 3, minor_breaks = NULL,
  labels = function(x) format(x, scientific = FALSE, big.mark = ",", trim = TRUE))
# ---- chord ----
p_chord <- ggchord(seq_data = seqs, ribbon_data = links, gene_data = genes) +
  geom_link_ribbon(fill = "#AFC7C0", alpha = .45) + geom_seq(show.legend = FALSE) + geom_gene(position = position_strand(.1), width = .04) +
  geom_seq_label(size = 4, seq_label_radius = .84, seq_label_hjust = .5) + coord_chord() +
  scale_seq_colour_manual(values = seq_palette) +
  scale_gene_fill_manual(values = gene_palette) + doc_theme
# ---- circular ----
data(plasmid_example_pUC19)
plasmid_seq <- plasmid_example_pUC19[c("accver", "length")]
feature_source <- find_common_features(plasmid_example_pUC19)
features <- feature_source[c("accver", "start", "end", "strand", "type", "anno")]

circle_features <- features[features$type %in% c("CDS", "rep_origin"), ]
p_circular <- ggchord(seq_data = plasmid_seq) + geom_seq() +
  geom_feature(aes(feature_fill = type), data = circle_features, position = position_plasmid()) +
  geom_feature_label_repel(data = circle_features, position = position_plasmid(), size = 3.3) +
  geom_seq_center_label() + scale_feature_fill_manual(values = feature_palette, breaks = unique(circle_features$type)) +
  coord_circular() + theme_ggchord_publication(axis.text = element_text(size = 9))
# ---- collinear ----
anchor_rows <- unlist(lapply(split(seq_len(nrow(genes)), genes$accver), function(i)
  i[which.max(genes$end[i] - genes$start[i])]))
anchors <- genes[anchor_rows, c("accver", "start", "end")]
p_collinear <- ggchord(seq_data = seqs, ribbon_data = links, gene_data = genes) +
  geom_link_ribbon(fill = "#AFC7C0", alpha = .4) + geom_seq() +
  geom_gene(position = position_strand(.1)) + geom_seq_label(size = 4) +
  coord_collinear(anchor = anchors, align = "center", orientation = "auto") + doc_axis() +
  scale_seq_colour_manual(values = seq_palette) +
  scale_gene_fill_manual(values = gene_palette) + doc_theme
# ---- adjustments ----
identity_links <- ribbon_data_example[c("qaccver", "saccver", "qstart", "qend", "sstart", "send", "pident")]
p_mapped <- ggchord(seq_data = seqs, ribbon_data = identity_links) + geom_seq() +
  geom_link_ribbon(aes(ribbon_fill = pident), alpha = .6) +
  scale_ribbon_fill_gradientn(colours = c("#E1EBE5", "#83A99B", "#315F59")) +
  scale_seq_colour_manual(values = seq_palette) +
  coord_chord(rotation = 20) + labs(title = "Sequence comparison") + doc_theme
# ---- labels ----
label_data <- transform(gene_data_example, strand = source_strand)
label_data <- label_data[label_data$anno != "hypothetical protein",
  c("accver", "start", "end", "strand", "anno")]
p_labels <- ggchord(seq_data = seqs, gene_data = genes) + geom_seq() + geom_gene() +
  geom_gene_label_repel(data = label_data, label_layout = "auto", size = 3.4, label_wrap = 22) +
  geom_seq_label(size = 4) + coord_collinear() + doc_axis() +
  scale_seq_colour_manual(values = seq_palette) + scale_gene_fill_manual(values = gene_palette) + doc_theme
# ---- tracks ----
data(plasmid_example_pUC19)
track_seq <- plasmid_example_pUC19[c("accver", "length")]
starts <- seq.int(1L, track_seq$length, by = 100L)
bins <- data.frame(accver = track_seq$accver, start = starts,
  end = pmin(starts + 99L, track_seq$length))
bins$gc <- vapply(seq_len(nrow(bins)), function(i) {
  bases <- strsplit(substr(plasmid_example_pUC19$sequence, bins$start[i], bins$end[i]), "")[[1]]
  mean(bases %in% c("G", "C"))
}, numeric(1))
p_tracks <- ggchord(seq_data = track_seq) + geom_seq() + geom_seq_center_label() +
  geom_track_bar(aes(track_value = gc, track_fill = gc), data = bins,
    position = position_track("gc", width = unit(6, "mm"))) +
  scale_track_value_continuous("gc", limits = c(0, 1), unit = "GC fraction") +
  scale_track_fill_gradientn("gc", colours = c("#EDF1EA", "#8AAE9C", "#345F58")) +
  coord_circular() + doc_theme
# ---- recognition ----
data(plasmid_example_pUC19)
target <- setNames(plasmid_example_pUC19$sequence, plasmid_example_pUC19$accver)
primers <- find_primer_bindings(target, set = "universal")
sites <- find_restriction_sites(target, enzymes = c("EcoRI", "BamHI", "HindIII"), methylation = "none")
p_sites <- ggchord(seq_data = plasmid_example_pUC19[c("accver", "length")]) + geom_seq() +
  geom_primer(data = primers, colour = "#7C83A2") +
  geom_restriction_site(data = sites, colour = "#455651") + coord_circular()
# ---- homologs ----
data(plasmid_example_pUC19)
data(plasmid_example_pETDuet_1)
homolog_source <- rbind(plasmid_example_pUC19, plasmid_example_pETDuet_1)
homolog_seq <- homolog_source[c("accver", "length")]
cluster_features <- do.call(rbind, lapply(seq_len(nrow(homolog_source)), function(i) {
  f <- find_common_features(homolog_source[i, ])
  f[f$type == "CDS" & f$anno == "AmpR", c("accver", "start", "end", "strand", "anno")]
}))
cluster_features$feature_id <- with(cluster_features, paste(accver, start, end, sep = ":"))
# Membership is computed from exact strand-corrected CDS DNA equality, not gene names.
cds_dna <- vapply(seq_len(nrow(cluster_features)), function(i) {
  f <- cluster_features[i, ]; dna <- homolog_source$sequence[match(f$accver, homolog_source$accver)]
  dna <- substr(dna, f$start, f$end)
  if (f$strand == "-") dna <- chartr("ACGT", "TGCA", paste(rev(strsplit(dna, "")[[1]]), collapse = ""))
  dna
}, "")
membership <- data.frame(feature_id = cluster_features$feature_id,
  cluster_id = paste0("exact_CDS_", match(cds_dna, unique(cds_dna))))
homologs <- link_ggchord_clusters(cluster_features, membership,
  pairs = data.frame(qaccver = homolog_seq$accver[1], saccver = homolog_seq$accver[2]))
loci <- find_ggchord_loci(homolog_seq, cluster_features,
  targets = cluster_features$feature_id[cluster_features$anno == "AmpR"][1], upstream = 100, downstream = 100)
p_homologs <- ggchord(seq_data = homolog_seq) +
  geom_link_ribbon(data = homologs$data, fill = "#ADC5B9", alpha = .6) + geom_seq() +
  geom_gene(data = cluster_features) + geom_seq_label(size = 4) +
  coord_collinear() + doc_axis() + scale_gene_fill_manual(values = gene_palette) + doc_theme
# ---- projection ----
# Run the exact-CDS comparison first; these are complete source-backed CDS.
coding <- transform(cluster_features, phase = 0)
cds <- build_ggchord_cds(coding, genetic_code = 11)
# This selects residues for a display example; it does not claim a detected domain.
domains <- data.frame(feature_id = "selected_residues", parent_id = coding$feature_id[1], start = 2, end = 20)
genomic <- project_ggchord_features(domains, cds)
# ---- focus ----
focus_loci <- data.frame(accver = seqs$accver,
  start = floor(seqs$length * .25), end = floor(seqs$length * .75))
focused <- focus_ggchord_data(seqs, links, genes, loci = focus_loci)
p_focus <- ggchord(seq_data = focused$seq_data, ribbon_data = focused$ribbon_data,
  gene_data = focused$gene_data) + geom_seq() +
  geom_link_ribbon(fill = "#ADC5B9", alpha = .5) + geom_gene() + geom_seq_break() +
  geom_seq_label(size = 4) + scale_seq_colour_manual(values = seq_palette) +
  scale_gene_fill_manual(values = gene_palette) +
  doc_axis("source") + coord_collinear() + doc_axis() + doc_theme
# ---- segment ----
windows <- data.frame(accver = rep(seqs$accver[1], 2),
  start = c(1, floor(seqs$length[1] * .65)), end = c(floor(seqs$length[1] * .4), seqs$length[1]),
  locus_id = c("left_window", "right_window"), row_id = seqs$accver[1], orientation = c(1, -1))
segmented <- segment_ggchord_data(seqs[1, ], windows = windows)
p_segments <- ggchord(seq_data = segmented$seq_data) + geom_seq() +
  coord_collinear(segments = segmented$segments, segment_gap = 2000)
# ---- validation ----
validation <- validate_ggchord_data(seqs, links, genes)
cleaned <- clean_ggchord_data(seqs, links, genes)
prepared <- prepare_ggchord_data(seqs, links, features = list(genes = genes))
# ---- provenance ----
geometry <- export_ggchord_layout(p_collinear, include = "seq")
rows <- unique(geometry$seq[c("layer_id", "component", "source_row")])
# ---- migration ----
p_migrated <- ggchord(plasmid_seq) + geom_seq() +
  geom_feature(data = features, shape = "block", width = .05) +
  geom_feature_label(data = features, label_rotation = 0) + coord_circular()
p_offsets <- ggchord(seqs, gene_data = genes) + geom_seq() +
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
sequence_source <- read_fasta_lengths(files = "examples/fasta/*.fna", source_file = TRUE)
seq_real <- sequence_source[c("accver", "length")]
annotation_source <- read_gff3(files = "examples/gff3/*.gff3", unstranded = "keep", source_file = TRUE)
feature_real <- annotation_source[c("accver", "start", "end", "strand", "type", "anno")]
alignment_source <- read_blast(files = "examples/blastn-standard/*.tsv", source_file = TRUE)
link_real <- alignment_source[c("qaccver", "saccver", "qstart", "qend", "sstart", "send", "pident", "length")]
real_validation <- validate_ggchord_data(seq_real, link_real)
# ---- realprepare ----
real_prepared <- prepare_ggchord_data(seq_real, link_real,
  features = list(annotation = feature_real))
real_filtered <- filter_ggchord_ribbons(link_real, min_pident = 80, min_length = 500)
# Explicit display subset: the six longest CDS per sequence; raw annotations stay intact.
display_rows <- unlist(lapply(split(seq_len(nrow(feature_real)), feature_real$accver),
  function(i) head(i[order(-(feature_real$end[i] - feature_real$start[i] + 1))], 6)))
real_display_features <- feature_real[display_rows, , drop = FALSE]
label_rows <- unlist(lapply(split(seq_len(nrow(real_display_features)), real_display_features$accver),
  function(i) head(i[real_display_features$anno[i] != "hypothetical protein"], 2)))
real_labels <- real_display_features[label_rows, , drop = FALSE]
# ---- realchord ----
p_real_chord <- ggchord(seq_data = seq_real, ribbon_data = real_filtered$data, validate = "none") +
  geom_link_ribbon(aes(ribbon_fill = pident), alpha = .65) + geom_seq(show.legend = FALSE) +
  geom_feature(aes(feature_fill = strand), data = real_display_features,
    width = .045, position = position_strand(.08)) +
  geom_feature_label_repel(data = real_labels, label_layout = "radial",
    label_wrap = 20, size = 3.4, position = position_strand(.08)) +
  geom_seq_label(size = 4, seq_label_radius = .84, seq_label_hjust = .5) + coord_chord() +
  scale_seq_colour_manual(values = seq_palette) +
  scale_feature_fill_manual(values = gene_palette, breaks = c("+", "-")) +
  scale_ribbon_fill_gradientn(colours = c("#DFEAE5", "#8BAA9C", "#395D55")) +
  doc_theme
# ---- reallinear ----
p_real_linear <- ggchord(seq_real, real_filtered$data, validate = "none") +
  geom_link_ribbon(fill = "#ADC5B9", alpha = .45) + geom_seq() + geom_seq_label(size = 4) +
  geom_feature(data = feature_real, feature_fill = "#71988A", show.legend = FALSE, position = position_strand(.18)) +
  geom_seq_scale(length = 10000) + doc_axis() +
  coord_collinear() + doc_axis() +
  scale_seq_colour_manual(values = seq_palette) + doc_theme
# ---- plasmid ----
data(plasmid_example_pUC19)
plasmid_seq <- plasmid_example_pUC19[c("accver", "length")]
plasmid_feature_source <- find_common_features(plasmid_example_pUC19)
# Dedicated primer layers already show bindings; avoid duplicate primer annotations.
plasmid_features <- plasmid_feature_source[plasmid_feature_source$type != "primer_bind",
  c("accver", "start", "end", "strand", "type", "anno")]
plasmid_primers <- find_primer_bindings(plasmid_example_pUC19, set = "universal")
plasmid_sites <- find_restriction_sites(plasmid_example_pUC19)
display_sites <- filter_restriction_sites(plasmid_sites, enzymes = c("EcoRI", "BamHI", "HindIII"))
primer_names <- primer_catalog(set = "universal")
plasmid_position <- position_feature_stack(spacing = .085, base_position = position_plasmid())
p_plasmid <- ggchord(seq_data = plasmid_seq) + geom_seq(seq_style = "double") +
  geom_feature_plasmid(aes(feature_fill = type), data = plasmid_features, position = plasmid_position) +
  geom_feature_label_repel(data = plasmid_features, position = plasmid_position,
    label_fit = "shrink", label_min_size = 2.4, size = 3.3) +
  geom_primer(data = plasmid_primers, colour = "#7C83A2") +
  geom_primer_label_repel(data = plasmid_primers, colour = "#656F92", max_overlaps = 0, size = 3) +
  geom_restriction_site(data = display_sites, colour = "#455651", label_size = 3) + geom_seq_center_label() +
  scale_feature_fill_manual(values = feature_palette) +
  coord_circular(rotation = 90) + theme_ggchord_publication(axis.text = element_text(size = 9))
# ---- cleaning ----
# Deliberately corrupt copies to demonstrate diagnostics; original rows remain intact.
problem_links <- ribbon_data_example[1:2, c("qaccver", "saccver", "qstart", "qend", "sstart", "send", "pident")]
problem_links$qaccver[2] <- "absent"
problem_links$pident[1] <- 103
problem_report <- validate_ggchord_data(seqs, problem_links)
repaired <- clean_ggchord_data(seqs, problem_links, unknown_id = "drop", invalid_pident = "clip")
# ---- hierarchy ----
# Each source record represents a separate complete phage genome.
contigs <- transform(seqs, genome = accver)
genomes <- data.frame(genome = seqs$accver)
hierarchy <- validate_ggchord_hierarchy(contigs, genomes, strict = TRUE)
# ---- adapters ----
renamed_sequences <- setNames(seqs, c("chrom", "size"))
sequence_input <- as_ggchord_seq(renamed_sequences, mapping = c(accver = "chrom", length = "size"))
renamed_features <- setNames(genes, c("chrom", "left", "right", "strand"))
feature_input <- as_ggchord_feature(renamed_features,
  mapping = c(accver = "chrom", start = "left", end = "right"))
link_input <- as_ggchord_link(links)
# ---- formats ----
# BED3 carries intervals only; the source annotations remain in features.
bed_path <- tempfile(fileext = ".bed")
bed_rows <- data.frame(accver = features$accver, start0 = features$start - 1L, end = features$end)
write.table(bed_rows, bed_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
bed_features <- read_bed(bed_path)
unlink(bed_path)
# Unchanged PAF rows from the public EMALE tutorial, not a simulated alignment.
paf_source <- read_paf("examples/format-inputs/emales-pair.paf")
paf_selected <- paf_source[!is.na(paf_source$mapq) & paf_source$mapq >= 30, ]
paf_links <- paf_selected[c("qaccver", "saccver", "qstart", "qend", "sstart", "send")]
paf_seqs <- unique(rbind(data.frame(accver = paf_source$qaccver, length = paf_source$qlen),
  data.frame(accver = paf_source$saccver, length = paf_source$slen)))
p_paf <- ggchord(seq_data = paf_seqs, ribbon_data = paf_links) +
  geom_link_ribbon(fill = "#ADC5BA", alpha = .5) + geom_seq(colour = "#547D73") +
  geom_seq_label(size = 4) + coord_collinear() + doc_axis() + doc_theme
# ---- standardblast ----
blast_source <- read_blast("examples/blastn-standard/MT108731.1__MT118296.1.tsv")
blast_links <- blast_source[c("qaccver", "saccver", "qstart", "qend", "sstart", "send")]
# A six-column table is sufficient for drawing; extra fields stay in blast_source.
custom_path <- tempfile(fileext = ".tsv")
write.table(blast_links, custom_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
custom_source <- read_blast(custom_path, format = "custom", col_names = names(blast_links))
custom_links <- custom_source[names(blast_links)]
unlink(custom_path)
# ---- protein ----
protein_source <- read_blast("examples/format-inputs/emales-mavirus-blastp.tsv", source_file = TRUE)
emale_source <- read_gff3("examples/format-inputs/emales-BVI_023A.gff3", source_file = TRUE)
# The source declares complete Prodigal CDS and translation table 11.
emale_coding <- emale_source[grepl("partial=00", emale_source$attributes),
  c("feature_id", "accver", "start", "end", "strand", "phase")]
emale_cds <- build_ggchord_cds(emale_coding, genetic_code = 11)
# Keep strong HSPs whose query protein has a known complete CDS.
# Gaps are retained: projecting one query span does not assert an ungapped pair.
protein_kept <- protein_source[protein_source$qaccver %in% emale_coding$feature_id &
  protein_source$evalue <= 1e-10 & protein_source$length >= 50, ]
protein_kept <- protein_kept[order(-protein_kept$bitscore), ]
selected_hsp <- protein_kept[1, ]
query_interval <- data.frame(feature_id = "query_HSP", parent_id = selected_hsp$qaccver,
  start = selected_hsp$qstart, end = selected_hsp$qend)
projected_query <- project_ggchord_features(query_interval, emale_cds)
# Subject projection additionally requires its own CDS map; it is not guessed from an accession.
# ---- genbank ----
# The original NCBI L09137.2 record identifies pUC19c; keep that identity.
genbank_result <- NULL
if (nzchar(Sys.which("python3")) && identical(suppressWarnings(system2("python3",
  c("-c", shQuote("from Bio import SeqIO")), stdout = FALSE, stderr = FALSE)), 0L)) {
  genbank_result <- read_genbank("examples/format-inputs/pUC19c-L09137.2.gb")
}
# ---- optionalranges ----
range_inputs <- NULL
if (requireNamespace("GenomicRanges", quietly = TRUE)) {
  unknown <- features[features$strand == ".", ][1, ]
  ranges <- GenomicRanges::GRanges(seqnames = unknown$accver,
    ranges = IRanges::IRanges(unknown$start, unknown$end), strand = "*",
    seqlengths = setNames(plasmid_seq$length, plasmid_seq$accver))
  range_inputs <- list(seq_data = as_ggchord_seq(ranges), features = as_ggchord_feature(ranges))
}
# ---- ribbons ----
unique_links <- deduplicate_ggchord_ribbons(real_filtered$data, by = "overlap",
  min_reciprocal_overlap = .9)
merged_links <- merge_ggchord_ribbons(unique_links$data, max_gap = 100,
  min_pident_difference = 2)
bundled_links <- bundle_ggchord_ribbons(link_real, seq_real, bins = 30)
p_density <- ggchord(seq_real, link_real, validate = "none") + geom_seq() +
  stat_ribbon_density(aes(ribbon_alpha = after_stat(density)), fill = "#71958A") +
  scale_ribbon_alpha_continuous(name = "Relative density", range = c(.1, .6),
    guide = guide_ggchord_legend(position = "bottom", nrow = 1)) +
  scale_seq_colour_manual(values = seq_palette) + coord_chord() + doc_theme
p_bundle <- ggchord(seq_real, link_real, validate = "none") + geom_seq() +
  stat_ribbon_bundle(bins = 30, fill = "#71958A") + coord_chord() + doc_theme
# ---- pointlinks ----
# Midpoints summarize existing HSPs for display; they are not new alignments.
points <- ribbon_data_example[c("qaccver", "saccver", "qstart", "qend", "sstart", "send", "pident")]
points$qpos <- floor((points$qstart + points$qend) / 2)
points$spos <- floor((points$sstart + points$send) / 2)
points <- points[c("qaccver", "saccver", "qpos", "spos", "pident")]
p_points <- ggchord(seq_data = seqs) + geom_seq() + geom_seq_label(size = 4) +
  geom_link_line(aes(link_colour = pident), data = points,
    arrow = arrow(length = unit(2, "mm")), linewidth = .7) +
  scale_link_colour_gradientn(colours = c("#CCDCD3", "#355F58"),
    name = "Identity (%)", guide = guide_ggchord_colourbar(available_aes = "link_colour",
      position = "bottom", direction = "horizontal",
      theme = theme(legend.key.width = unit(35, "mm"), legend.key.height = unit(3, "mm")))) +
  scale_seq_colour_manual(values = seq_palette) + coord_chord() + doc_theme
# ---- featurestack ----
stack_features <- features[c("accver", "start", "end", "strand", "type", "anno")]
stack_position <- position_feature_stack(spacing = .12)
feature_shapes <- setNames(ifelse(names(feature_palette) == "CDS", "arrow", "block"), names(feature_palette))
p_stack <- ggchord(seq_data = plasmid_seq) + geom_seq() +
  geom_feature(aes(feature_fill = type, feature_shape = type), data = stack_features,
    position = stack_position) +
  geom_feature_label_repel(data = stack_features, position = stack_position, size = 3.2) +
  scale_feature_shape_manual(values = feature_shapes) +
  scale_feature_fill_manual(values = feature_palette) + coord_collinear() + doc_axis() + doc_theme
# ---- labellayouts ----
# Compare layouts with at most the first and last informative label per sequence.
comparison_labels <- do.call(rbind, lapply(split(label_data, label_data$accver),
  function(x) x[unique(c(1L, nrow(x))), , drop = FALSE]))
label_plots <- lapply(c("auto", "radial", "arc"), function(layout)
  ggchord(seq_data = seqs, gene_data = genes) + geom_seq() + geom_gene() +
    geom_gene_label_repel(data = comparison_labels, label_layout = layout, size = 3.4, label_wrap = 22) +
    geom_seq_label(size = 4) + scale_seq_colour_manual(values = seq_palette) +
    scale_gene_fill_manual(values = gene_palette) +
    coord_chord(xlim = c(-1.8, 1.8), ylim = c(-1.6, 1.6)) +
    theme_ggchord_publication(axis = element_blank()) + labs(title = layout))
names(label_plots) <- c("auto", "radial", "arc")
p_manual_labels <- ggchord(seq_data = seqs, gene_data = genes) + geom_seq() + geom_gene() +
  geom_gene_label(data = label_data, label_rotation = 30) + coord_collinear()
# ---- decorations ----
regions <- data.frame(accver = links$qaccver[1], start = min(links$qstart[1], links$qend[1]),
  end = max(links$qstart[1], links$qend[1]), category = "alignment interval")
p_decorated <- ggchord(seq_data = seqs, ribbon_data = links) +
  geom_seq_region(aes(region_fill = category), data = regions, colour = NA) +
  geom_link_ribbon(fill = "#ADC5B9", alpha = .4) + geom_seq() + geom_seq_label(size = 4) +
  scale_region_fill_manual(values = c("alignment interval" = "#E3D0A5")) +
  scale_seq_colour_manual(values = seq_palette) + geom_seq_scale(length = 10000) +
  coord_collinear() + doc_axis() + doc_theme
# ---- multitrack ----
line_bins <- transform(bins, pos = floor((start + end) / 2),
  gc_mean = as.numeric(stats::filter(gc, rep(1 / 3, 3), sides = 2)))
line_bins <- line_bins[is.finite(line_bins$gc_mean), c("accver", "pos", "gc_mean")]
categories <- features[c("accver", "start", "end", "type")]
p_multitrack <- ggchord(seq_data = track_seq) + geom_seq() + geom_seq_center_label() +
  geom_track_bar(aes(track_value = gc, track_fill = gc), data = bins,
    position = position_track("gc", width = unit(6, "mm"))) +
  geom_track_line(aes(track_value = gc_mean), data = line_bins,
    position = position_track("gc", width = unit(6, "mm")), colour = "#454D65") +
  scale_track_value_continuous("gc", limits = c(0, 1), unit = "GC fraction") +
  scale_track_fill_gradientn("gc", colours = c("#EDF1EA", "#8AAE9C", "#345F58")) +
  geom_track_tile(aes(track_fill = type), data = categories,
    position = position_track("features", order = 2, width = unit(2, "mm"))) +
  scale_track_fill_manual("features", values = feature_palette) + coord_circular() + doc_theme
focused_bins <- focus_ggchord_track(bins, track_seq,
  loci = data.frame(accver = track_seq$accver, start = 500, end = 1500), type = "interval")
# ---- linkprojection ----
# Exact CDS DNA equality supplies evidence for these selected ungapped residues.
matching <- homologs$data[1, ]
paired_coding <- transform(cluster_features, phase = 0)
paired_cds <- build_ggchord_cds(paired_coding)
protein_link <- data.frame(qfeature_id = matching$qfeature_id, sfeature_id = matching$sfeature_id,
  qstart = 1, qend = 20, sstart = 1, send = 20)
projected_links <- project_ggchord_links(protein_link, paired_cds, ungapped = TRUE)
# ---- mappingposition ----
source_position <- focus_loci$start[1] + 100
source_mark <- map_ggchord_position(p_focus, seqs$accver[1], source_position, source = TRUE)
p_marked <- p_focus + ggplot2::geom_point(data = source_mark, aes(x, y),
  inherit.aes = FALSE, colour = "#A16E7D", size = 3)
# ---- styling ----
p_styled <- ggchord(seq_data = seqs, ribbon_data = links, gene_data = genes) + geom_seq() +
  geom_link_ribbon(aes(ribbon_fill = qaccver), alpha = .5) + geom_gene() +
  geom_seq_label(size = 4) + coord_collinear() + doc_axis() +
  scale_seq_colour_manual(values = seq_palette) + scale_gene_fill_manual(values = gene_palette) +
  scale_ribbon_fill_manual(values = setNames(c("#D4E3DC", "#E5DCCB", "#D8DFE9", "#E4D8DC"), seqs$accver)) +
  guides(gene_fill = guide_ggchord_legend(nrow = 1)) + doc_theme
theme_variants <- list(minimal = theme_ggchord_minimal(), dark = theme_ggchord_dark(), plasmid = theme_ggchord_plasmid())
# ---- layoutresults ----
full_layout <- get_chord_layout(p_collinear)
layer_tables <- export_ggchord_layout(p_collinear, include = c("seq", "gene"), original_data = TRUE)
layout_choice <- optimize_ggchord_layout(seqs, links, optimize = c("order", "orientation"), weight = "count")
# ---- parameters ----
# Sequence IDs remain stable when the data is reordered; positional lists do not.
radii <- setNames(c(1.1, 1, 1, 1), seqs$accver)
widths <- c("+" = .045, "-" = .065)
p_parameters <- ggchord(seq_data = seqs, gene_data = genes) +
  geom_seq(seq_radius = radii) + geom_gene(width = widths,
    position = position_strand(c("+" = .08, "-" = -.1))) +
  geom_seq_label(size = 4) + scale_seq_colour_manual(values = seq_palette) +
  scale_gene_fill_manual(values = gene_palette) + coord_chord() + doc_theme
# ---- end ----
