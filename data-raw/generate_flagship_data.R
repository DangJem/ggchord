# Rebuild the ready-to-plot flagship datasets from unchanged bundled sources.
# Run after data-raw/generate_website_data.R, from the package root.
# No ggplot objects, device measurements or label placement solutions are saved.
pkgload::load_all(".", quiet = TRUE)
library(dplyr)
show_sequences <- genome_comparison_example$seq_data
show_links <- genome_comparison_example$links
show_cross <- genome_comparison_example$cross_links
chord_genes <- genome_comparison_example$chord_features
show_linear_features <- genome_comparison_example$linear_features
show_adjacent <- genome_comparison_example$adjacent_links
# Preserve the four display directions selected for the published Chord.
# Orientation belongs to the presentation; alignment endpoints remain in source order.
chord_orientation <- c("MT108731.1" = -1, "MT118296.1" = 1,
  "OQ646790.1" = 1, "OR222515.1" = -1)
radial_orientation <- setNames(rep(1, nrow(show_sequences)), show_sequences$accver)

radial_seqs <- seq_data_example |>
  select(accver, length)
radial_hsp <- ribbon_data_example |>
  mutate(source_hsp_id = paste0("bundled:", row_number())) |>
  filter(length >= 900, pident >= 75) |>
  group_by(qaccver, saccver) |>
  slice_max(length, n = 2, with_ties = FALSE) |>
  ungroup()
radial_features <- gene_data_example |>
  transmute(accver, start, end, strand = source_strand, anno, type = "CDS")
radial_events <- radial_hsp |>
  transmute(accver = qaccver, pos = qstart)
flagship_comparison_example <- list(
  comparison = list(seq_data = radial_seqs, links = radial_hsp, features = radial_features, events = radial_events),
  chord = list(seq_data = show_sequences, links = show_links,
    features = chord_genes, labels = chord_genes,
    orientation = chord_orientation),
  collinear = list(seq_data = show_sequences, adjacent_links = show_adjacent,
    cross_links = show_cross, features = show_linear_features,
    labels = arrange(show_linear_features, desc(nchar(anno)), accver, start)),
  radial = list(seq_data = show_sequences, links = show_links,
    features = chord_genes, labels = chord_genes, orientation = radial_orientation),
  sankey = list(seq_data = show_sequences, adjacent_links = show_adjacent,
    cross_links = show_cross, features = show_linear_features, labels = show_linear_features)
)
flagship_sequences <- genome_comparison_example$seq_data
flagship_hsp <- genome_comparison_example$source_alignments
flagship_features <- filter(genome_comparison_example$annotations, type == "CDS")
dotplot_pair_counts <- flagship_hsp |>
  count(qaccver, saccver, name = "HSP")
dotplot_gap <- 4000
dotplot_frames <- flagship_sequences |>
  mutate(offset = c(0, head(cumsum(length + dotplot_gap), -1L)),
    min = offset + 1, max = offset + length)
dotplot_cells <- tidyr::expand_grid(
  qaccver = flagship_sequences$accver, saccver = flagship_sequences$accver
) |>
  left_join(dotplot_pair_counts, by = c("qaccver", "saccver")) |>
  mutate(HSP = coalesce(HSP, 0L)) |>
  left_join(dotplot_frames |>
    transmute(qaccver = accver, x = min, xmax = max), by = "qaccver") |>
  left_join(dotplot_frames |>
    transmute(saccver = accver, y = max, ymin = min), by = "saccver") |>
  mutate(x = x + (xmax - x) * 0.04, y = y - (y - ymin) * 0.05,
    label = if_else(HSP > 0L, paste(HSP, "HSP"), "No input HSP"))

flagship_dotplot_example <- list(seq_data = flagship_sequences,
  links = flagship_hsp, cell_labels = dotplot_cells, pair_counts = dotplot_pair_counts,
  axis_gap = dotplot_gap)
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
flagship_focused <- focus_ggchord_data(
  ggchord_data(flagship_sequences, features = list(genes = flagship_features),
    links = list(observed = flagship_selected)),
  mutate(flagship_windows, locus_id = paste0(accver, "_detail")))
fusion_detail_features <- flagship_focused$features$genes |>
  group_by(accver) |>
  arrange(anno == "hypothetical protein", start, .by_group = TRUE) |>
  slice_head(n = 2) |>
  ungroup()
fusion_context_sequences <- flagship_sequences |> filter(accver %in% flagship_windows$accver)
fusion_detail_source_features <- flagship_features |>
  semi_join(fusion_detail_features, by = join_by(accver == source_accver, feature_id))
# Select the second displayed CDS on E's lower sequence, then read its real bases.
fusion_base_cds <- fusion_detail_source_features |>
  filter(accver == tail(flagship_focused$sequences$source_accver, 1)) |>
  arrange(start) |>
  slice(2)
fusion_base_genome <- genome_comparison_example$sequences |>
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

flagship_fusion_example <- list(
  seq_data = flagship_sequences, links = fusion_display_hsp,
  features = fusion_display_features, orientation = fusion_orientation,
  all_features = flagship_features,
  events = transmute(flagship_hsp, accver = qaccver, pos = pmin(qstart, qend)),
  relation_ab = mutate(fusion_relation_ab, relation = "A - B"),
  relation_ac = mutate(fusion_relation_ac, relation = "A - C"),
  relation_cd = mutate(flagship_selected, relation = "C - D"),
  selected = flagship_selected, windows = flagship_windows, focused = flagship_focused,
  detail_features = fusion_detail_features, context_sequences = fusion_context_sequences,
  context_features = fusion_context_features, base_cds = fusion_base_cds,
  base_start = fusion_base_start, base_end = fusion_base_end,
  base_cells = fusion_base_cells, focus_regions = fusion_focus_regions
)

# Dense pETDuet-1 flagship: keep its existing unique-site parent collection.
sequence <- mutate(plasmid_example_pETDuet_1, topology = "circular")
dense_sites <- find_restriction_sites(sequence)
flagship_circular_example <- list(dense = list(
  seq_data = select(sequence, accver, length),
  features = find_common_features(sequence, mode = "reference"),
  primers = find_primer_bindings(sequence, set = "reference"),
  all_site_count = nrow(dense_sites),
  sites = filter_restriction_sites(dense_sites, min_site_length = 6, cuts = 1,
    commercial_only = TRUE)
))
# Complete reference maps and the older three-cut gallery are distinct selections.
objects <- c("plasmid_example_pSB1C3", "plasmid_example_pUC19", "plasmid_example_pBR322",
  "plasmid_example_pBluescript_II_SK_plus", "plasmid_example_pET_28a_plus",
  "plasmid_example_pETDuet_1", "plasmid_example_pDONR221", "plasmid_example_pcDNA3_1_plus",
  "plasmid_example_pTRE_Tight_BI", "plasmid_example_pTRIPZ", "plasmid_example_pCAMBIA1300",
  "plasmid_example_pEarleyGate_201", "plasmid_example_pSpCas9_BB_2A_GFP_PX458")
records <- lapply(objects, function(object) {
  e <- new.env(parent = emptyenv())
  data(list = object, package = "ggchord", envir = e)
  sequence <- mutate(e[[object]], topology = "circular")
  features <- find_common_features(sequence, mode = "reference")
  all_sites <- find_restriction_sites(sequence)
  profile <- attr(all_sites, "reference_enzyme_profiles")
  # The catalogue reproduces the reference display selection. Group equivalent
  # recognition starts by preferred enzyme, retaining custom reference sets.
  sites <- ggchord_filter_restriction_profile(all_sites, set = "reference",
    parent_set = "commercial_nonredundant")
  policy <- "Saved reference profile; original catalogue display representatives"
  if (nrow(profile) && profile$profile_type == "none") {
    sites <- ggchord_filter_restriction_profile(all_sites, set = "unique_6plus",
      parent_set = "commercial_nonredundant")
    policy <- "Explicit unique 6+ profile; saved None does not reproduce the reference PNG"
  }
  targets <- sequence$length * c(1, 3, 5) / 6
  nearest <- vapply(targets, function(x) which.min(abs(all_sites$pos - x)), integer(1))
  list(seq_data = select(sequence, accver, length), source_sequence = sequence,
    features = features, primers = find_primer_bindings(sequence, set = "reference"),
    all_site_count = nrow(all_sites), sites = sites, policy = policy, gallery_sites = all_sites[unique(nearest), ])
})
names(records) <- vapply(records, function(x) x$seq_data$accver, character(1))
flagship_circular_example$references <- records
flagship_circular_example$summary <- bind_rows(lapply(records, function(x)
  tibble(accver = x$seq_data$accver, length = x$seq_data$length, features = nrow(x$features))))
# Store plain source/display tables and explicit presets, never compiled plots.
for (name in c("flagship_comparison_example", "flagship_dotplot_example",
               "flagship_fusion_example", "flagship_circular_example")) {
  destination <- file.path("data", paste0(name, ".rda"))
  temporary <- tempfile(fileext = ".rda")
  save(list = name, file = temporary, compress = "xz")
  stopifnot(file.copy(temporary, destination, overwrite = TRUE))
  unlink(temporary)
}
inputs <- c("data/genome_comparison_example.rda", "R/sysdata.rda",
  paste0("data/", objects, ".rda"))
write.table(data.frame(file = inputs, md5 = unname(tools::md5sum(inputs))),
  "data-raw/flagship-data-sources.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
