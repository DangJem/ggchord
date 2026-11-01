# Shared executable documentation snippets, v0.16.5.
# Marker names are consumed by tools/build-website.py and tools/sync-tutorials.py.
# All A/B/C data and signal values below are synthetic teaching inputs.
# ---- setup ----
library(ggchord)
seqs <- data.frame(accver = c("A", "B", "C"), length = c(1000, 800, 1200))
links <- data.frame(qaccver = c("A", "B"), saccver = c("B", "C"),
  qstart = c(100, 200), qend = c(350, 500),
  sstart = c(150, 400), send = c(400, 700), pident = c(96, 89))
genes <- data.frame(accver = c("A", "B", "C"),
  start = c(150, 200, 400), end = c(300, 350, 650),
  strand = c("+", "+", "-"), anno = c("target", "target", "target"))
# ---- chord ----
p_chord <- ggchord(seqs, links, genes) +
  geom_seq() + geom_link_ribbon() + geom_gene() + geom_seq_label() +
  coord_chord()
# ---- circular ----
features <- data.frame(accver = "A", start = c(80, 380, 700),
  end = c(280, 610, 850), strand = c("+", "-", "."),
  type = c("CDS", "CDS", "region"), label = c("alpha", "beta", "origin"))
p_circular <- ggchord(seqs[1, ]) + geom_seq() +
  geom_feature(data = features) +
  geom_feature_label_repel(data = features) + coord_circular() +
  theme_ggchord(axis = element_blank())
# ---- collinear ----
p_collinear <- ggchord(seqs, links, genes) +
  geom_link_ribbon(alpha = .3) + geom_seq() +
  geom_gene(position = position_strand(.1)) + geom_seq_label() +
  coord_collinear(anchor = genes[c("accver", "start", "end")],
    align = "center", orientation = "auto")
# ---- adjustments ----
p_mapped <- ggchord(seqs, links) + geom_seq() +
  geom_link_ribbon(aes(ribbon_fill = pident)) +
  scale_ribbon_fill_gradientn(colours = c("#E7EFEA", "#246D68")) +
  coord_chord(rotation = 20) + labs(title = "Sequence comparison") +
  theme_ggchord() + theme(legend.position = "bottom")
# ---- labels ----
p_labels <- ggchord(seqs, gene_data = genes) + geom_seq() +
  geom_gene() + geom_gene_label_repel(label_layout = "auto", size = 3) +
  coord_collinear()
# ---- tracks ----
bins <- data.frame(accver = "A", start = seq(1, 951, 50),
  end = seq(50, 1000, 50), signal = sin(seq(0, 2 * pi, length.out = 20)))
p_tracks <- ggchord(seqs[1, ]) + geom_seq() +
  geom_track_bar(aes(track_value = signal, track_fill = signal), data = bins,
    position = position_track("signal", width = unit(5, "mm"))) +
  scale_track_value_continuous("signal", limits = c(-1, 1), unit = "arbitrary units") +
  scale_track_fill_gradientn("signal", colours = c("#31537D", "#F3F4EE", "#D26B38")) +
  coord_circular()
# ---- recognition ----
target <- c(circle = "TTGCAAAAAGAATTCACG")
primers <- find_primer_bindings(target, primers = c(test = "ACGTTGCA"),
  min_annealed_bases = 8)
sites <- find_restriction_sites(target, enzymes = "EcoRI", methylation = "none")
p_sites <- ggchord(data.frame(accver = "circle", length = nchar(target))) +
  geom_seq() + geom_primer(data = primers) +
  geom_restriction_site(data = sites) + coord_circular()
# ---- projection ----
coding <- data.frame(feature_id = "cdsA", accver = "A", start = c(100, 200),
  end = c(103, 204), strand = "+", phase = c(0, 2))
cds <- build_ggchord_cds(coding, genetic_code = 11)
domains <- data.frame(feature_id = "domainA", parent_id = "cdsA", start = 2, end = 3)
genomic <- project_ggchord_features(domains, cds)
# ---- focus ----
focused <- focus_ggchord_data(seqs, links, genes,
  loci = data.frame(accver = c("A", "B", "C"),
    start = c(50, 100, 300), end = c(600, 600, 950)))
p_focus <- ggchord(focused$seq_data, focused$ribbon_data, focused$gene_data) +
  geom_seq() + geom_link_ribbon() + geom_gene() + geom_seq_break() +
  scale_seq_position_continuous(coordinates = "source") + coord_collinear()
# ---- segment ----
windows <- data.frame(accver = c("A", "A"), start = c(50, 650),
  end = c(400, 950), locus_id = c("A_left", "A_right"),
  row_id = "A", orientation = c(1, -1))
segmented <- segment_ggchord_data(seqs[1, ], windows = windows)
p_segments <- ggchord(segmented$seq_data) + geom_seq() +
  coord_collinear(segments = segmented$segments, segment_gap = 80)
# ---- validation ----
validation <- validate_ggchord_data(seqs, links, genes)
cleaned <- clean_ggchord_data(seqs, links, genes)
prepared <- prepare_ggchord_data(seqs, links, features = list(genes = genes))
# ---- provenance ----
geometry <- export_ggchord_layout(p_collinear, include = "seq")
rows <- unique(geometry$seq[c("layer_id", "component", "source_row")])
# ---- migration ----
p_migrated <- ggchord(seqs[1, ]) + geom_seq() +
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
link_files <- read_blast(files = "examples/blastn/*.o7", source_file = TRUE)
feature_files <- read_gff3(files = "examples/gff3/*.gff3", unstranded = "keep", source_file = TRUE)
# ---- realimport ----
seq_real <- read_fasta_lengths(files = "examples/fasta/*.fna", source_file = TRUE)
feature_real <- read_gff3(files = "examples/gff3/*.gff3", unstranded = "keep", source_file = TRUE)
link_real <- read_blast(files = "examples/blastn/*.o7", source_file = TRUE)
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
  function(i) head(i, 2)))
real_labels <- real_display_features[label_rows, , drop = FALSE]
# ---- realchord ----
p_real_chord <- ggchord(seq_real, real_filtered$data, validate = "none") +
  geom_link_ribbon(aes(ribbon_fill = pident), alpha = .35) + geom_seq() +
  geom_feature(aes(feature_fill = strand), data = real_display_features,
    width = .045, position = position_strand(.08)) +
  geom_feature_label_repel(data = real_labels, label_layout = "radial",
    label_wrap = 20, size = 2.5, position = position_strand(.08)) +
  geom_seq_label() + coord_chord() + theme_ggchord_publication() +
  theme(legend.position = "bottom")
# ---- reallinear ----
p_real_linear <- ggchord(seq_real, real_filtered$data, validate = "none") +
  geom_link_ribbon(alpha = .25) + geom_seq() + geom_seq_label() +
  geom_feature(data = feature_real, position = position_strand(.08)) +
  geom_seq_scale(length = 10000) + coord_collinear() +
  theme_ggchord_publication() + theme(legend.position = "bottom")
# ---- plasmid ----
data(plasmid_example_pUC19)
plasmid_features <- find_common_features(plasmid_example_pUC19)
plasmid_primers <- find_primer_bindings(plasmid_example_pUC19, set = "universal")
plasmid_sites <- find_restriction_sites(plasmid_example_pUC19)
display_sites <- filter_restriction_sites(plasmid_sites,
  enzymes = c("EcoRI", "BamHI", "HindIII"))
primer_names <- primer_catalog(set = "universal")
plasmid_position <- position_feature_stack(spacing = .085,
  base_position = position_plasmid())
p_plasmid <- ggchord(plasmid_example_pUC19) + geom_seq(seq_style = "double") +
  geom_feature_plasmid(data = plasmid_features, position = plasmid_position) +
  geom_feature_label_repel(data = plasmid_features, position = plasmid_position,
    label_fit = "shrink", label_min_size = 1.8) +
  geom_primer(data = plasmid_primers) +
  geom_primer_label_repel(data = plasmid_primers, max_overlaps = 0) +
  geom_restriction_site(data = display_sites) + geom_seq_center_label() +
  coord_circular(rotation = 90) + theme_ggchord_plasmid()
# ---- cleaning ----
problem_links <- data.frame(qaccver = c("A", "absent"), saccver = "B",
  qstart = c(80, 20), qend = c(20, 40), sstart = 50, send = 120, pident = c(103, 90))
problem_report <- validate_ggchord_data(seqs, problem_links)
repaired <- clean_ggchord_data(seqs, problem_links, unknown_id = "drop",
  reversed_interval = "sort", invalid_pident = "clip")
# Deliberate teaching repairs, not a recommendation to clip scientific values silently.
# ---- hierarchy ----
contigs <- transform(seqs, genome = c("sample1", "sample1", "sample2"))
genomes <- data.frame(genome = c("sample1", "sample2"))
hierarchy <- validate_ggchord_hierarchy(contigs, genomes, strict = TRUE)
# ---- adapters ----
sequence_input <- as_ggchord_seq(data.frame(chrom = "A", size = 1000),
  mapping = c(accver = "chrom", length = "size"))
feature_input <- as_ggchord_feature(data.frame(chrom = "A", left = 10,
  right = 100, strand = "+"), mapping = c(accver = "chrom", start = "left", end = "right"))
link_input <- as_ggchord_link(links)
# ---- formats ----
bed_path <- tempfile(fileext = ".bed")
writeLines("A\t0\t100\tdemo\t0\t-", bed_path)
bed_features <- read_bed(bed_path)
paf_path <- tempfile(fileext = ".paf")
writeLines("A\t1000\t0\t100\t-\tB\t800\t200\t300\t95\t100\t255", paf_path)
paf_links <- read_paf(paf_path)
unlink(c(bed_path, paf_path))
# ---- genbank ----
# Optional reader: the selected Python must already contain Biopython.
genbank_result <- NULL
if (nzchar(Sys.which("python3")) && identical(suppressWarnings(system2("python3",
  c("-c", shQuote("from Bio import SeqIO")), stdout = FALSE, stderr = FALSE)), 0L)) {
  gb_path <- tempfile(fileext = ".gb")
  writeLines(c("LOCUS       demo                      12 bp    DNA     linear   BCT 01-JAN-2000",
    "DEFINITION  Synthetic teaching record.", "ACCESSION   demo", "VERSION     demo.1",
    "FEATURES             Location/Qualifiers", "     misc_feature    2..8",
    "                     /note=\"site\"", "ORIGIN", "        1 acgtacgtacgt", "//"), gb_path)
  genbank_result <- read_genbank(gb_path)
  unlink(gb_path)
}
# ---- optionalranges ----
range_inputs <- NULL
if (requireNamespace("GenomicRanges", quietly = TRUE)) {
  ranges <- GenomicRanges::GRanges(seqnames = "A", ranges = IRanges::IRanges(10, 100),
    strand = "*", seqlengths = c(A = 1000L))
  range_inputs <- list(seq_data = as_ggchord_seq(ranges), features = as_ggchord_feature(ranges))
}
# ---- ribbons ----
unique_links <- deduplicate_ggchord_ribbons(real_filtered$data, by = "overlap",
  min_reciprocal_overlap = .9)
merged_links <- merge_ggchord_ribbons(unique_links$data, max_gap = 100,
  min_pident_difference = 2)
bundled_links <- bundle_ggchord_ribbons(link_real, seq_real, bins = 30)
p_density <- ggchord(seq_real, link_real, validate = "none") + geom_seq() +
  stat_ribbon_density(aes(ribbon_alpha = after_stat(density))) + coord_chord()
p_bundle <- ggchord(seq_real, link_real, validate = "none") + geom_seq() +
  stat_ribbon_bundle(bins = 30) + coord_chord()
# ---- pointlinks ----
points <- data.frame(qaccver = c("A", "A"), saccver = c("B", "C"),
  qpos = c(200, 200), spos = c(400, 600), support = c(4, 10))
p_points <- ggchord(seqs) + geom_seq() + geom_seq_label() +
  geom_link_line(aes(link_colour = support, link_linewidth = support), data = points,
    arrow = arrow(length = unit(2, "mm"))) +
  scale_link_colour_gradientn(colours = c("#BBCFC8", "#265E55"), guide = guide_ggchord_colourbar(available_aes = "link_colour")) +
  scale_link_linewidth_continuous(range = c(.3, 1)) + coord_chord()
# ---- featurestack ----
stack_features <- data.frame(feature_id = paste0("f", 1:4), accver = "A",
  start = c(100, 150, 300, 650), end = c(350, 450, 550, 800),
  strand = c("+", "+", "-", "."), type = c("CDS", "region", "CDS", "region"),
  label = c("coding alpha", "overlap", "coding beta", "unstranded"))
stack_position <- position_feature_stack(spacing = .12)
p_stack <- ggchord(seqs[1, ]) + geom_seq() +
  geom_feature(aes(feature_fill = type, feature_shape = type), data = stack_features,
    position = stack_position) +
  geom_feature_label_repel(data = stack_features, position = stack_position) +
  scale_feature_shape_manual(values = c(CDS = "arrow", region = "block")) +
  scale_feature_fill_manual(values = c(CDS = "#4C8C7C", region = "#D3A05A")) + coord_collinear()
# ---- labellayouts ----
label_plots <- lapply(c("auto", "radial", "arc"), function(layout)
  ggchord(seqs, gene_data = genes) + geom_seq() + geom_gene() +
    geom_gene_label_repel(label_layout = layout, size = 3) + coord_chord() +
    labs(title = layout))
names(label_plots) <- c("auto", "radial", "arc")
p_manual_labels <- ggchord(seqs, gene_data = genes) + geom_seq() + geom_gene() +
  geom_gene_label(label_rotation = 30) + coord_collinear()
# ---- decorations ----
regions <- data.frame(accver = "A", start = 80, end = 420, category = "selected")
p_decorated <- ggchord(seqs, links) +
  geom_seq_region(aes(region_fill = category), data = regions) +
  geom_link_ribbon(alpha = .25) + geom_seq() + geom_seq_label() +
  scale_region_fill_manual(values = c(selected = "#E1B566")) +
  geom_seq_scale(length = 200) + coord_collinear()
# ---- multitrack ----
line_bins <- transform(bins, pos = floor((start + end) / 2), signal = signal / 2)
categories <- data.frame(accver = "A", start = c(100, 500, 800),
  end = c(250, 700, 950), category = c("one", "two", "one"))
p_multitrack <- ggchord(seqs[1, ]) + geom_seq() +
  geom_track_bar(aes(track_value = signal, track_fill = signal), data = bins,
    position = position_track("signal", width = unit(5, "mm"))) +
  geom_track_line(aes(track_value = signal), data = line_bins,
    position = position_track("signal", width = unit(5, "mm")), colour = "#244B68") +
  scale_track_value_continuous("signal", limits = c(-1, 1), unit = "arbitrary units") +
  scale_track_fill_gradientn("signal", colours = c("#31537D", "#F3F4EE", "#D26B38")) +
  geom_track_tile(aes(track_fill = category), data = categories,
    position = position_track("sites", side = "outside", order = 2, width = unit(2, "mm"))) +
  scale_track_fill_manual("sites", values = c(one = "#6E9E88", two = "#DCBA66")) +
  coord_circular()
focused_bins <- focus_ggchord_track(bins, seqs[1, ],
  loci = data.frame(accver = "A", start = 100, end = 700), type = "interval")
# ---- homologs ----
cluster_features <- transform(genes, feature_id = paste0("gene", 1:3))
membership <- data.frame(feature_id = cluster_features$feature_id, cluster_id = "demo_cluster")
homologs <- link_ggchord_clusters(cluster_features, membership,
  pairs = data.frame(qaccver = c("A", "B"), saccver = c("B", "C")))
loci <- find_ggchord_loci(seqs, cluster_features, targets = "gene1", upstream = 50, downstream = 100)
p_homologs <- ggchord(seqs) + geom_link_ribbon(data = homologs$data, fill = "#9CBBAF") +
  geom_seq() + geom_gene(data = genes) + geom_seq_label() + coord_collinear()
# ---- linkprojection ----
paired_coding <- data.frame(feature_id = c("a", "b"), accver = c("A", "B"),
  start = c(10, 100), end = c(18, 108), strand = c("+", "-"), phase = 0)
paired_cds <- build_ggchord_cds(paired_coding)
protein_link <- data.frame(qfeature_id = "a", sfeature_id = "b",
  qstart = 1, qend = 2, sstart = 2, send = 3)
projected_links <- project_ggchord_links(protein_link, paired_cds, ungapped = TRUE)
# ---- mappingposition ----
source_mark <- map_ggchord_position(p_focus, "A", 200, source = TRUE)
p_marked <- p_focus + ggplot2::geom_point(data = source_mark, aes(x, y),
  inherit.aes = FALSE, colour = "#A23B40", size = 2)
# ---- styling ----
p_styled <- ggchord(seqs, links, genes) + geom_seq() +
  geom_link_ribbon(aes(ribbon_fill = qaccver), alpha = .3) + geom_gene() +
  geom_seq_label() + coord_collinear() +
  scale_seq_colour_manual(values = c(A = "#487A6A", B = "#B68048", C = "#637DA0")) +
  scale_gene_fill_manual(values = c("+" = "#B68048", "-" = "#637DA0")) +
  scale_ribbon_fill_manual(values = c("A" = "#D4E5D8", "B" = "#E3D9C8", "C" = "#D4DFEB")) +
  guides(gene_fill = guide_ggchord_legend(nrow = 1)) +
  theme_ggchord_publication() + theme(legend.position = "bottom")
theme_variants <- list(minimal = theme_ggchord_minimal(), dark = theme_ggchord_dark(),
  plasmid = theme_ggchord_plasmid())
# ---- layoutresults ----
full_layout <- get_chord_layout(p_collinear)
layer_tables <- export_ggchord_layout(p_collinear, include = c("seq", "gene"), original_data = TRUE)
layout_choice <- optimize_ggchord_layout(seqs, links, optimize = c("order", "orientation"), weight = "count")
# ---- end ----
