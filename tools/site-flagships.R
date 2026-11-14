# Shared source-backed website figures; source this file after loading ggchord.
if (!exists("p_show_chord") || !exists("p_show_collinear")) {
  shared_lines <- readLines("tools/documentation-examples.R", warn = FALSE)
  markers <- grep("^# ---- [a-z]+ ----$", shared_lines)
  shared_blocks <- setNames(lapply(seq_len(length(markers) - 1L), function(i)
    shared_lines[seq.int(markers[i] + 1L, markers[i + 1L] - 1L)]),
    sub("^# ---- ([a-z]+) ----$", "\\1", shared_lines[markers[-length(markers)]]))
  for (key in c("sequencepalette", "docstyle", "showcasedata", "chordshowcasegenes", "showcasechord", "showcasecollinear"))
    invisible(eval(parse(text = shared_blocks[[key]])))
}
# ---- flagshipdata ----
library(ggchord)
library(dplyr)
flagship_paper <- "#FFFCF5"
flagship_sequences <- read_fasta_lengths("examples/fasta/*.fna") |>
  select(accver, length)
flagship_hsp <- read_blast(files = "examples/blastn-standard/*.tsv", source_file = TRUE) |>
  mutate(source_hsp_id = paste0("source-HSP-", row_number()))
flagship_features <- read_gff3("examples/gff3/*.gff3", unstranded = "keep") |>
  filter(type == "CDS")
flagship_events <- flagship_hsp |>
  transmute(accver = qaccver, pos = pmin(qstart, qend))
flagship_colours <- c(
  "MT108731.1" = "#E64B35", "MT118296.1" = "#233E4A",
  "OQ646790.1" = "#269B82", "OR222515.1" = "#D59B17"
)
flagship_theme <- theme_ggchord_publication() +
  theme(
    plot.background = element_rect(fill = flagship_paper, colour = NA),
    panel.background = element_rect(fill = flagship_paper, colour = NA),
    legend.background = element_rect(fill = flagship_paper, colour = NA),
    legend.key = element_rect(fill = flagship_paper, colour = NA),
    legend.position = "right", legend.box = "vertical", legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    plot.margin = margin(12, 16, 12, 16)
  )
# ---- flagshipmulti ----
relationship_identity <- ggplot2::guide_coloursteps(
  available_aes = "ribbon_fill", even.steps = FALSE, show.limits = TRUE,
  position = "bottom", direction = "horizontal",
  theme = theme(
    legend.key.width = unit(130, "mm"), legend.key.height = unit(4, "mm"),
    legend.title = element_text(size = 10), legend.text = element_text(size = 9),
    legend.margin = margin(2, 2, 2, 2)
  )
)
relationship_guides <- guides(
  seq_colour = "none",
  gene_fill = guide_ggchord_legend(
    nrow = 1, position = "bottom", override.aes = list(strand = c("+", "-"))
  ),
  ribbon_fill = relationship_identity
)
relationship_theme <- doc_theme + theme(
  legend.position = "bottom", legend.box = "vertical",
  legend.title = element_text(size = 10, face = "bold"),
  legend.text = element_text(size = 9),
  legend.box.spacing = unit(3, "mm"), legend.spacing.y = unit(2, "mm"),
  legend.margin = margin(2, 2, 2, 2),
  ggchord.axis.text = element_text(size = 9),
  plot.margin = margin(8, 10, 8, 10),
  plot.title = element_text(size = 17), plot.subtitle = element_text(size = 10)
)
radial_labels <- chord_genes |>
  group_by(accver) |>
  arrange(grepl("hypothetical", anno, ignore.case = TRUE), desc(end - start), .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup()
radial_position <- position_strand(c("+" = -0.055, "-" = 0.055))
flag_radial <- ggchord(show_sequences, show_links) +
  geom_link_ribbon(adjacent_only = FALSE, link_avoid = "smooth", alpha = 0.46, linewidth = 0.20) +
  geom_seq(seq_orientation = show_optimization$seq_orientation, seq_gap = 0.3, linewidth = 1.1) +
  geom_gene(data = chord_genes, position = radial_position, width = 0.065) +
  geom_gene_label_repel(data = radial_labels, label_layout = "radial", label_wrap = 22,
    size = 3.3, label_max_lines = 3, position = radial_position) +
  geom_seq_label(size = 4.4, seq_label_orientation = "horizontal") +
  scale_seq_colour_manual(values = seq_palette, name = "Sequence") +
  scale_gene_fill_manual(values = gene_palette, breaks = c("+", "-"), name = "Strand") +
  scale_seq_position_continuous(n.breaks = 2, minor_breaks = NULL,
    labels = function(x) format(x, scientific = FALSE, big.mark = ",", trim = TRUE)) +
  coord_radial_synteny(start_angle = 45, inner_radius = 0.7) +
  relationship_theme + relationship_guides +
  labs(title = "Four genomes, a radial view", subtitle = "Original Chord sources; one representative CDS name per sequence")
sankey_labels <- show_linear_features |>
  group_by(accver) |>
  arrange(grepl("hypothetical", anno, ignore.case = TRUE), desc(end - start), .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup()
flag_sankey <- ggchord(show_sequences, gene_data = show_linear_features) +
  geom_link_ribbon(aes(ribbon_fill = pident), data = show_adjacent,
    link_avoid = "smooth", alpha = 0.48, linewidth = 0.25) +
  geom_link_ribbon(aes(ribbon_fill = pident), data = show_cross,
    adjacent_only = FALSE, link_avoid = "smooth", alpha = 0.28, linewidth = 0.3) +
  geom_seq(seq_gap = 0.30, linewidth = 1.1) +
  geom_gene(position = position_strand(0.09), width = 0.075) +
  geom_gene_label(data = sankey_labels,
    label_orientation = "horizontal", label_side = "outside", label_overlap = "allow",
    label_radial_offset = 0.25, label_wrap = 16, position = position_strand(0.09), size = 3.3) +
  geom_seq_label(size = 4.4) +
  scale_seq_colour_manual(values = seq_palette, name = "Sequence") +
  scale_gene_fill_manual(values = gene_palette, breaks = c("+", "-"), name = "Strand") +
  scale_seq_position_continuous(n.breaks = 2, minor_breaks = NULL,
    labels = function(x) format(x, scientific = FALSE, big.mark = ",", trim = TRUE)) +
  coord_sankey(column_gap = 5.2) + relationship_theme + relationship_guides +
  labs(title = "Four genomes, mirrored relationships", subtitle = "Original Collinear sources; independent query and subject displays")
# ---- fusionoverview ----
flagship_position <- position_strand(0.045)
flagship_base <- ggchord(flagship_sequences, flagship_hsp, feature_data = flagship_features) +
  geom_link_ribbon(adjacent_only = FALSE, alpha = 0.45, linewidth = 0.12, allow_repeated = TRUE) +
  geom_seq(aes(seq_colour = accver), linewidth = 1.1) +
  geom_gene(aes(gene_fill = strand), position = flagship_position, width = 0.025) +
  geom_seq_label(size = 3.2) +
  scale_seq_colour_manual(values = flagship_colours, name = "Sequence") +
  scale_gene_fill_manual(values = c("+" = "#CE1536", "-" = "#346B9D", "." = "#8291A6"), name = "Biological strand") +
  flagship_theme
fusion_chord <- flagship_base + coord_chord()
fusion_collinear <- flagship_base + coord_collinear()
fusion_radial <- flagship_base + coord_radial_synteny()
# ---- flagshipdot ----
flag_dotplot <- ggchord(flagship_sequences, flagship_hsp) +
  geom_alignment_segment(aes(link_colour = pident, link_linewidth = length), alpha = 0.85) +
  scale_link_linewidth_continuous(range = c(0.45, 1.7), name = "HSP length (bp)") +
  coord_dotplot(x_seqs = flagship_sequences$accver, y_seqs = flagship_sequences$accver,
    gap = 4000) +
  labs(x = "Query position (bp)", y = "Subject position (bp)") +
  flagship_theme +
  theme(panel.grid.major = element_line(colour = "#DFDED5", linewidth = 0.25))
# ---- flagshipfusion ----
flagship_selected <- flagship_hsp |> slice(1)
flagship_windows <- bind_rows(
  flagship_selected |> transmute(accver = qaccver, start = pmin(qstart, qend), end = pmax(qstart, qend)),
  flagship_selected |> transmute(accver = saccver, start = pmin(sstart, send), end = pmax(sstart, send))
)
flagship_focused <- focus_ggchord_data(flagship_sequences, ribbon_data = flagship_selected,
  feature_data = flagship_features, loci = flagship_windows)
flagship_detail <- ggchord(flagship_focused$seq_data, flagship_focused$ribbon_data,
  feature_data = flagship_focused$feature_data) +
  geom_link_ribbon(allow_repeated = TRUE, alpha = 0.45) +
  geom_seq(aes(seq_colour = accver)) +
  geom_gene(aes(gene_fill = strand), width = 0.045) +
  geom_seq_label(size = 2.6) + coord_collinear() +
  scale_seq_colour_manual(values = flagship_colours) +
  scale_gene_fill_manual(values = c("+" = "#CE1536", "-" = "#346B9D", "." = "#8291A6")) +
  flagship_theme + theme(legend.position = "none")
flagship_quiet <- theme(legend.position = "none")
flag_fusion <- ggchord_fusion(fusion_chord + flagship_quiet,
  center_box = c(0.02, 0.46, 0.53, 0.98)) +
  fusion_panel(fusion_collinear + flagship_quiet, id = "linear", box = c(0.57, 0.98, 0.53, 0.98)) +
  fusion_panel(fusion_radial + flagship_quiet, id = "radial", box = c(0.02, 0.46, 0.02, 0.43)) +
  fusion_panel(flagship_detail, id = "detail", box = c(0.57, 0.98, 0.02, 0.43)) +
  fusion_ribbon(flagship_hsp |> slice(2:5), qpanel = "center", spanel = "linear",
    fill = "#269B82", alpha = 0.28, allow_repeated = TRUE) +
  fusion_ribbon(flagship_hsp |> slice(6:9), qpanel = "center", spanel = "radial",
    fill = "#D59B17", alpha = 0.28, allow_repeated = TRUE) +
  fusion_ribbon(flagship_selected, qpanel = "center", spanel = "detail",
    fill = "#E64B35", alpha = 0.32, allow_repeated = TRUE) +
  theme(plot.background = element_rect(fill = flagship_paper, colour = NA))
# ---- referencebuilder ----
make_reference_map <- function(sequence, paper = "white") {
  features <- find_common_features(sequence)
  all_sites <- find_restriction_sites(sequence)
  profile <- attr(all_sites, "reference_enzyme_profiles")
  sites <- filter_restriction_sites(all_sites, set = "reference",
    parent_set = "commercial_nonredundant")
  policy <- "Saved reference profile; commercial nonredundant representatives"
  if (nrow(profile) && profile$profile_type == "none") {
    sites <- filter_restriction_sites(all_sites, set = "unique_6plus",
      parent_set = "commercial_nonredundant")
    policy <- "Explicit unique 6+ profile; saved None does not reproduce the reference PNG"
  }
  primers <- find_primer_bindings(sequence, set = "reference")
  position <- position_feature_stack(spacing = 0.085, base_position = position_plasmid())
  plot <- ggchord(sequence[c("accver", "length")], feature_data = features) +
    geom_seq(seq_style = "double") +
    geom_feature_plasmid(position = position) +
    geom_feature_label_repel(position = position, label_fit = "shrink", label_min_size = 2.5) +
    geom_seq_center_label() + coord_circular() + theme_ggchord_plasmid()
  if (nrow(sites)) plot <- plot + geom_restriction_site(data = sites, label_size = 2.5)
  if (nrow(primers)) plot <- plot + geom_primer(data = primers) +
    geom_primer_label_repel(data = primers, size = 2.5)
  if (paper != "white") plot <- plot + theme(
    plot.background = element_rect(fill = paper, colour = NA),
    panel.background = element_rect(fill = paper, colour = NA))
  list(plot = plot, features = features, primers = primers, sites = sites,
    all_sites = all_sites, policy = policy)
}
# ---- flagshipcircle ----
flagship_circular_source <- make_reference_map(plasmid_example_pETDuet_1, paper = flagship_paper)
flag_circular <- flagship_circular_source$plot
# ---- flagend ----
