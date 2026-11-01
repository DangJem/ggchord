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
# ---- end ----
