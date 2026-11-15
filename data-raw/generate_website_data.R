# Rebuild the complete genome comparison and sequence-detail teaching data.
# Run Rscript data-raw/generate_website_data.R from the package root.
# Source records and coordinates are unchanged. Selection rules below reproduce
# the existing website flagships, including row order and biological strands.
pkgload::load_all(".", quiet = TRUE)
library(dplyr)
# Unchanged repository inputs; the full tables and display subsets are saved separately.
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

genome_comparison_example <- list(
  sequences = read_fasta_lengths("examples/fasta/*.fna", keep_sequence = TRUE),
  seq_data = show_sequences,
  annotations = show_annotations,
  alignments = show_alignments,
  source_alignments = read_blast("examples/blastn-standard/*.tsv", source_file = TRUE) |>
    mutate(source_hsp_id = paste0("source-HSP-", row_number())),
  genes = show_genes, links = show_links, features = show_features,
  chord_features = chord_genes, linear_features = show_linear_features,
  adjacent_links = show_adjacent, cross_links = show_cross
)
sequence_detail_example <- list(
  plasmid_source = read_fasta_lengths("examples/plasmid/pUC19.fna", keep_sequence = TRUE, source_file = TRUE),
  paf = read_paf("examples/format-inputs/emales-pair.paf"),
  dna_alignment = read_aligned_fasta("examples/alignment/opuntia-aligned.fasta", "DNA", source_file = TRUE),
  protein_alignment = read_aligned_fasta("examples/alignment/protein-aligned.fasta", "AA"),
  aligned_strings = read_fasta_lengths("examples/alignment/opuntia-aligned.fasta", keep_sequence = TRUE),
  protein_hits = read_blast("examples/format-inputs/emales-mavirus-blastp.tsv", source_file = TRUE),
  source_coding_features = read_gff3("examples/format-inputs/emales-BVI_023A.gff3", unstranded = "keep", source_file = TRUE),
  coding_genome = read_fasta_lengths("examples/format-inputs/emales-BVI_023A.fna", keep_sequence = TRUE),
  coding_features = read_gff3("examples/format-inputs/emales-BVI_023A.gff3", unstranded = "keep"),
  coding_proteins = read_fasta_lengths("examples/format-inputs/protein-pair.faa", keep_sequence = TRUE)
)
# Import provenance retains source file names and coordinate conventions.
# A separate manifest fingerprints the unchanged source files for rebuilding.
inputs <- c("examples/format-inputs/emales-mavirus-blastp.tsv", "examples/plasmid/pUC19.fna", Sys.glob("examples/fasta/*.fna"), Sys.glob("examples/gff3/*.gff3"),
  Sys.glob("examples/blastn-standard/*.tsv"),
  "examples/format-inputs/emales-pair.paf", "examples/alignment/opuntia-aligned.fasta",
  "examples/alignment/protein-aligned.fasta", "examples/format-inputs/emales-BVI_023A.fna",
  "examples/format-inputs/emales-BVI_023A.gff3", "examples/format-inputs/protein-pair.faa")
write.table(data.frame(file = inputs, md5 = unname(tools::md5sum(inputs))),
  "data-raw/website-data-sources.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
save(genome_comparison_example, file = "data/genome_comparison_example.rda", compress = "xz")
save(sequence_detail_example, file = "data/sequence_detail_example.rda", compress = "xz")
