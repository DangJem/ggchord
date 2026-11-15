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
library(ggplot2)
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
  seq_colour = guide_ggchord_legend(nrow = 1, position = "bottom"),
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
radial_labels <- chord_genes
radial_orientation <- show_optimization$seq_orientation
radial_orientation[c("MT108731.1", "OR222515.1")] <- 1
# Adjustable annular lift for non-neighbours, in layout units; longest spoke = 2.
radial_cross_bulge <- 0.20
radial_position <- position_strand(c("+" = -0.055, "-" = 0.055))
flag_radial <- ggchord(show_sequences, show_links) +
  geom_link_ribbon(adjacent_only = FALSE, link_avoid = "smooth", alpha = 0.46, linewidth = 0.20) +
  geom_seq(seq_orientation = radial_orientation, seq_gap = 0.3, linewidth = 1.1) +
  geom_gene(data = chord_genes, position = radial_position, width = 0.065) +
  geom_gene_label_repel(data = radial_labels, label_side = "auto", label_wrap = 20,
    label_max_lines = 3, position = radial_position, size = 3.1) +
  geom_seq_label(size = 4.4, seq_label_orientation = "horizontal",
    seq_label_vjust = c("MT108731.1" = 0.5, "MT118296.1" = 0.5,
      "OQ646790.1" = 1.5, "OR222515.1" = 0.5)) +
  scale_seq_colour_manual(values = seq_palette, name = "Sequence") +
  scale_gene_fill_manual(values = gene_palette, breaks = c("+", "-"), name = "Strand") +
  scale_seq_position_continuous(n.breaks = 2, minor_breaks = NULL,
    labels = function(x) format(x, scientific = FALSE, big.mark = ",", trim = TRUE)) +
  coord_radial_synteny(start_angle = 45, inner_radius = 0.7,
    nonadjacent_bulge = radial_cross_bulge) +
  relationship_theme + relationship_guides +
  labs(title = "Four genomes, a radial view", subtitle = "Original Chord sources; all 20 CDS annotations")
sankey_labels <- show_linear_features
flag_sankey <- ggchord(show_sequences, gene_data = show_linear_features) +
  geom_link_ribbon(aes(ribbon_fill = pident), data = show_adjacent,
    link_avoid = "smooth", alpha = 0.48, linewidth = 0.25) +
  geom_link_ribbon(aes(ribbon_fill = pident), data = show_cross,
    adjacent_only = FALSE, link_avoid = "smooth", alpha = 0.28, linewidth = 0.3) +
  geom_seq(seq_gap = 0.30, linewidth = 1.1) +
  geom_gene(position = position_strand(0.09), width = 0.075) +
  geom_gene_label_repel(data = sankey_labels, label_wrap = 20,
    label_max_lines = 3, position = position_strand(0.09), size = 3.1) +
  geom_seq_label(size = 4.4) +
  scale_seq_colour_manual(values = seq_palette, name = "Sequence") +
  scale_gene_fill_manual(values = gene_palette, breaks = c("+", "-"), name = "Strand") +
  scale_seq_position_continuous(n.breaks = 2, minor_breaks = NULL,
    labels = function(x) format(x, scientific = FALSE, big.mark = ",", trim = TRUE)) +
  coord_sankey(column_gap = 4.2, ratio = 1.15) + relationship_theme + relationship_guides +
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
dotplot_pair_counts <- flagship_hsp |>
  count(qaccver, saccver, name = "HSP")
dotplot_frames <- export_ggchord_layout(
  ggchord(flagship_sequences, flagship_hsp) + geom_alignment_segment() +
    coord_dotplot(x_seqs = flagship_sequences$accver,
      y_seqs = flagship_sequences$accver, gap = 4000)
)$axis_layout
dotplot_cells <- tidyr::expand_grid(
  qaccver = flagship_sequences$accver, saccver = flagship_sequences$accver
) |>
  left_join(dotplot_pair_counts, by = c("qaccver", "saccver")) |>
  mutate(HSP = coalesce(HSP, 0L)) |>
  left_join(dotplot_frames |> filter(axis == "x") |>
    transmute(qaccver = accver, x = min, xmax = max), by = "qaccver") |>
  left_join(dotplot_frames |> filter(axis == "y") |>
    transmute(saccver = accver, y = max, ymin = min), by = "saccver") |>
  mutate(x = x + (xmax - x) * 0.04, y = y - (y - ymin) * 0.05,
    label = if_else(HSP > 0L, paste(HSP, "HSP"), "No input HSP"))
flag_dotplot <- ggchord(flagship_sequences, flagship_hsp) +
  geom_alignment_segment(aes(link_colour = pident, link_linewidth = length), alpha = 0.85) +
  scale_link_linewidth_continuous(range = c(0.8, 1.8), name = "HSP length (bp)") +
  coord_dotplot(x_seqs = flagship_sequences$accver, y_seqs = flagship_sequences$accver,
    gap = 4000) +
  ggplot2::geom_text(data = dotplot_cells, aes(x = x, y = y, label = label),
    inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.8, colour = "#566E67") +
  labs(title = "Six observed sequence pairs · 32 HSPs",
    subtitle = "One query→subject direction per pair; empty cells contain no input records") +
  flagship_theme +
  theme(
    panel.grid = element_blank(), legend.position = "bottom", legend.box = "horizontal",
    axis.text = element_text(size = 9),
    axis.text.x.top = element_text(size = 10, face = "bold", margin = margin(b = 9)),
    axis.text.y.left = element_text(size = 10, face = "bold", margin = margin(r = 9)),
    axis.text.y.right = element_text(size = 9),
    axis.ticks.length = unit(1.5, "mm"),
    plot.margin = margin(12, 16, 12, 16)
  )
# ---- flagshipfusion ----
fusion_orientation <- setNames(rep(1, nrow(flagship_sequences)), flagship_sequences$accver)
fusion_display_hsp <- flagship_hsp |>
  group_by(qaccver, saccver) |>
  arrange(abs(length - 4000), .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup()
fusion_display_features <- flagship_features |>
  group_by(accver) |>
  arrange(start, .by_group = TRUE) |>
  mutate(display_bin = ntile(row_number(), 3)) |>
  group_by(accver, display_bin) |>
  arrange(anno == "hypothetical protein", desc(end - start), .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  select(-display_bin)
fusion_relation_ab <- fusion_display_hsp |> filter(qaccver == "MT108731.1", saccver == "MT118296.1")
fusion_relation_ac <- fusion_display_hsp |> filter(qaccver == "OQ646790.1", saccver == "OR222515.1")
flagship_selected <- fusion_display_hsp |> filter(qaccver == "MT118296.1", saccver == "OQ646790.1")
flagship_windows <- bind_rows(
  flagship_selected |> transmute(accver = qaccver, start = pmin(qstart, qend), end = pmax(qstart, qend)),
  flagship_selected |> transmute(accver = saccver, start = pmin(sstart, send), end = pmax(sstart, send))
) |>
  left_join(flagship_sequences, by = "accver") |>
  mutate(start = pmax(1, start - 800), end = pmin(length, end + 800)) |>
  select(accver, start, end)
flagship_focused <- focus_ggchord_data(flagship_sequences, ribbon_data = flagship_selected,
  feature_data = flagship_features, loci = flagship_windows)
fusion_detail_features <- flagship_focused$feature_data |>
  group_by(accver) |>
  arrange(anno == "hypothetical protein", start, .by_group = TRUE) |>
  slice_head(n = 2) |>
  ungroup()
fusion_strand_colours <- c("+" = "#CE1536", "-" = "#346B9D")
fusion_shared_styles <- list(
  scale_seq_colour_manual(values = flagship_colours, limits = names(flagship_colours), drop = FALSE, name = "Sequence"),
  scale_gene_fill_manual(values = fusion_strand_colours, limits = c("+", "-"), drop = FALSE, name = "Strand"),
  guides(seq_colour = guide_ggchord_legend(nrow = 2, position = "bottom"),
    gene_fill = guide_ggchord_legend(nrow = 1, position = "bottom", override.aes = list(strand = c("+", "-"))),
    ribbon_fill = ggplot2::guide_coloursteps(available_aes = "ribbon_fill", even.steps = FALSE,
      show.limits = TRUE, position = "bottom", direction = "horizontal",
      theme = theme(legend.key.width = grid::unit(85, "mm"), legend.key.height = grid::unit(3, "mm")))),
  flagship_theme + theme(legend.position = "bottom", legend.box = "horizontal",
    legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8),
    legend.margin = margin(2, 2, 2, 2), legend.spacing.x = grid::unit(2, "mm"),
    plot.margin = margin(3, 5, 3, 5), plot.title = element_text(size = 15, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9, colour = "#52646A"),
    plot.title.position = "plot", ggchord.axis.text = element_text(size = 9))
)
make_fusion_view <- function(coord, title, subtitle, label_radius = NULL) {
  label_layer <- if (is.null(label_radius)) {
    geom_seq_label(size = 3.5, seq_label_orientation = "horizontal")
  } else {
    geom_seq_label(size = 3.5, seq_label_radius = label_radius,
      seq_label_orientation = "horizontal")
  }
  ggchord(flagship_sequences, fusion_display_hsp, feature_data = fusion_display_features) +
    geom_link_ribbon(adjacent_only = FALSE, alpha = 0.25, linewidth = 0.10,
      nonadjacent_linewidth = 0.22, allow_repeated = TRUE) +
    geom_seq(aes(seq_colour = accver), seq_orientation = fusion_orientation, linewidth = 1.1) +
    geom_gene(aes(gene_fill = strand), position = flagship_position, width = 0.04) +
    label_layer +
    coord + fusion_shared_styles + labs(title = title, subtitle = subtitle)
}
fusion_chord_view <- make_fusion_view(coord_chord(rotation = 45,
    xlim = c(-1.22, 1.02), ylim = c(-1.02, 1.22), fit = "manual"),
  "A  |  Chord", "6 HSPs; 3 CDS per genome", label_radius = 1.12) +
  theme(ggchord.axis.text = element_blank(), ggchord.axis.ticks = element_blank())
fusion_linear_view <- make_fusion_view(coord_collinear(ratio = 0.78),
  "B  |  Collinear", "The same intervals in linear coordinates") +
  scale_seq_position_continuous(n.breaks = 2, minor_breaks = NULL)
fusion_radial_view <- make_fusion_view(coord_radial_synteny(start_angle = 45, inner_radius = 0.7),
  "C  |  Radial", "Shared intervals; A / D") +
  theme(ggchord.axis.text = element_blank(), ggchord.axis.ticks = element_blank())
fusion_zoom_view <- ggchord(flagship_focused$seq_data, flagship_focused$ribbon_data,
  feature_data = fusion_detail_features) +
  geom_link_ribbon(allow_repeated = TRUE, alpha = 0.25, linewidth = 0.10) +
  geom_seq(aes(seq_colour = accver), seq_orientation = 1, linewidth = 1.1) +
  geom_gene(aes(gene_fill = strand), position = flagship_position, width = 0.04) +
  geom_gene_label_repel(size = 3.05, label_wrap = 22, label_max_lines = 2, position = flagship_position) +
  geom_seq_label(size = 3.5, seq_label_orientation = "horizontal") +
  coord_collinear(ratio = 0.85, fit = "manual",
    xlim = c(-0.42, 2.08), ylim = c(-0.72, 0.22)) + fusion_shared_styles +
  scale_seq_position_continuous(n.breaks = 3, minor_breaks = NULL, coordinates = "source") +
  labs(title = "E  |  Local zoom", subtitle = "Source windows + 800 bp flanks") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
fusion_context_sequences <- flagship_sequences |> filter(accver %in% flagship_windows$accver)
fusion_detail_source_features <- flagship_features |>
  semi_join(fusion_detail_features, by = c("accver", "feature_id"))
# Select the second displayed CDS on E's lower sequence, then read its real bases.
fusion_base_cds <- fusion_detail_source_features |>
  filter(accver == tail(flagship_focused$seq_data$accver, 1)) |>
  arrange(start) |>
  slice(2)
fusion_base_genome <- read_fasta_lengths("examples/fasta/*.fna", keep_sequence = TRUE) |>
  filter(accver == fusion_base_cds$accver)
fusion_base_start <- if (fusion_base_cds$strand == "+") fusion_base_cds$start else fusion_base_cds$end - 59
fusion_base_end <- fusion_base_start + 59
fusion_base_string <- substring(fusion_base_genome$sequence, fusion_base_start, fusion_base_end)
if (fusion_base_cds$strand == "-") {
  fusion_base_string <- paste0(chartr("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN",
    rev(strsplit(fusion_base_string, "", fixed = TRUE)[[1]])), collapse = "")
}
fusion_base_cells <- as_ggchord_alignment(setNames(fusion_base_string, fusion_base_cds$feature_id),
  "DNA", source_data = data.frame(accver = fusion_base_cds$feature_id,
    source_accver = fusion_base_cds$accver,
    start = if (fusion_base_cds$strand == "+") fusion_base_start else fusion_base_end,
    strand = fusion_base_cds$strand, sequence = fusion_base_genome$sequence))
fusion_base_view <- ggplot(fusion_base_cells, aes(alignment_column, accver, label = symbol)) +
  geom_residue(overflow = "error", show.legend = FALSE) +
  scale_x_continuous(breaks = seq(1, 60, 10),
    labels = function(x) fusion_base_cells$source_position[match(x, fusion_base_cells$alignment_column)],
    expand = expansion(add = 0.6)) +
  scale_y_discrete(labels = setNames("DNA", fusion_base_cds$feature_id), expand = expansion(add = 0.7)) +
  labs(title = "F  |  Nucleotide detail", x = "Source position (bp)", y = NULL,
    subtitle = paste0(fusion_base_cds$accver, " / ", fusion_base_cds$feature_id,
      "  |  first 60 coding bases (", fusion_base_cds$strand, ")  |  ",
      fusion_base_start, "-", fusion_base_end, " bp")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), plot.background = element_rect(fill = flagship_paper, colour = NA),
    axis.text = element_text(colour = "#25363B"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, colour = "#52646A", hjust = 0.5), plot.margin = margin(4, 12, 4, 12))
fusion_base_view <- register_ggchord_alignment(fusion_base_view, fusion_base_cells)
fusion_focus_regions <- bind_rows(
  flagship_windows |> mutate(source_panel = "center", detail_panel = "zoom"),
  data.frame(accver = fusion_base_cds$accver, start = fusion_base_start, end = fusion_base_end,
    source_panel = "zoom", detail_panel = "bases")
)
fusion_context_features <- fusion_display_features |>
  inner_join(flagship_windows, by = "accver", suffix = c("", ".window")) |>
  filter(end < start.window | start > end.window) |>
  select(-start.window, -end.window) |>
  bind_rows(fusion_detail_source_features) |>
  distinct(accver, feature_id, .keep_all = TRUE)
fusion_detail_view <- ggchord(fusion_context_sequences, flagship_selected, feature_data = fusion_context_features) +
  geom_link_ribbon(allow_repeated = TRUE, alpha = 0.25, linewidth = 0.10) +
  geom_seq(aes(seq_colour = accver), seq_orientation = 1, linewidth = 1.1) +
  geom_gene(aes(gene_fill = strand), position = flagship_position, width = 0.04) +
  geom_seq_label(size = 3.5, seq_label_orientation = "horizontal") +
  coord_collinear(ratio = 0.7) + fusion_shared_styles +
  scale_seq_position_continuous(n.breaks = 2, minor_breaks = NULL) +
  labs(title = "D  |  Two-genome context", subtitle = "Full-length sequences; outlined windows are enlarged in E")
flagship_detail <- fusion_zoom_view
fusion_connections <- list(
  fusion_ribbon(fusion_relation_ab |> mutate(relation = "A - B"), aes(fill = relation),
    qpanel = "chord", spanel = "linear", alpha = 0.65, colour = "#3E5763", linewidth = 0.35,
    curvature = 0.09, allow_repeated = TRUE),
  fusion_ribbon(fusion_relation_ac |> mutate(relation = "A - C"), aes(fill = relation),
    qpanel = "chord", spanel = "radial", alpha = 0.65, colour = "#3E5763", linewidth = 0.35,
    curvature = 0.2, allow_repeated = TRUE),
  fusion_ribbon(flagship_selected |> mutate(relation = "C - D"), aes(fill = relation),
    qpanel = "radial", spanel = "center", alpha = 0.65, colour = "#3E5763", linewidth = 0.35,
    curvature = 0.2, route = "avoid", route_width = 1.5, allow_repeated = TRUE)
)
make_flagship_fusion <- function(preview = FALSE) {
  output_styles <- list(
    guides(ribbon_fill = ggplot2::guide_coloursteps(available_aes = "ribbon_fill",
      even.steps = FALSE, show.limits = TRUE, position = "bottom", direction = "horizontal",
      theme = theme(legend.key.width = grid::unit(if (preview) 90 else 85, "mm"),
        legend.key.height = grid::unit(3, "mm")))),
    theme(plot.title = element_text(size = if (preview) 15 else 16),
      plot.subtitle = element_text(size = if (preview) 9 else 10),
      legend.text = element_text(size = if (preview) 8 else 9))
  )
  lift <- function(box) c(box[1:2], 0.20 + 0.80 * box[3:4])
  out <- magnify_ggchord(fusion_detail_view + output_styles, region = fusion_focus_regions,
    detail_plot = list(zoom = fusion_zoom_view + output_styles, bases = fusion_base_view), guides = "collect",
    center_box = lift(c(0.465, 0.97, 0.385, 0.59)),
    detail_box = list(zoom = lift(c(0.49, 0.99, 0.01, 0.36)), bases = c(0.025, 0.985, 0, 0.17)),
    colour = "#725094", alpha = 0.20, highlight_width = 2.2,
    linewidth = 0.40, connector = "envelope", route = "avoid") +
    fusion_panel(fusion_chord_view + output_styles, id = "chord",
      box = lift(c(0.01, 0.455, if (preview) 0.52 else 0.54, 0.99))) +
    fusion_panel(fusion_linear_view + output_styles, id = "linear", box = lift(c(0.52, 0.995, 0.60, 0.99))) +
    fusion_panel(fusion_radial_view + output_styles, id = "radial",
      box = lift(c(0.01, 0.455, 0, if (preview) 0.50 else 0.52)))
  for (connection in fusion_connections) out <- out + connection
  out <- out + scale_fill_manual(name = if (preview) "Links" else "Correspondence", values = c("A - B" = "#527FB8",
    "A - C" = "#CE9741", "C - D" = "#9064AC")) +
    guides(fill = guide_legend(nrow = 1, position = "bottom",
      theme = if(preview) theme(legend.key.width = grid::unit(3, "mm"),
        legend.key.height = grid::unit(3, "mm")) else NULL)) +
    theme(plot.background = element_rect(fill = flagship_paper, colour = NA),
      legend.position = "bottom", legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8), plot.title = element_text(size = 22, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 11, colour = "#52646A"),
      plot.caption = element_text(size = 10, colour = "#52646A", hjust = 0),
      plot.margin = margin(8, 10, 8, 10))
  if (!preview) out <- out + labs(title = "One dataset, six linked views",
    subtitle = "Six selected HSPs and twelve source CDS per overview; full inputs retained separately",
    caption = "Coloured ribbons: observed HSPs. Source-window leaders connect D to E; the selected CDS bases connect E to F.")
  out
}
flag_fusion <- make_flagship_fusion()
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
