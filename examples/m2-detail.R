# Run from the source checkout. Optional dependencies: dplyr and Biostrings.
# Source: examples/format-inputs/README.md, unchanged EMALE genome/GFF3/protein.
library(ggchord)
library(dplyr)
genome <- read_fasta_lengths(
  'examples/format-inputs/emales-BVI_023A.fna', keep_sequence = TRUE
)
coding <- read_gff3('examples/format-inputs/emales-BVI_023A.gff3') |>
  filter(feature_id == 'BVI_023A_0009')
# Source annotation is complete, negative strand, phase 0, start_type=ATG.
# Table 11 is an explicit bacterial-code choice, checked against the source protein.
cds <- build_ggchord_cds(coding, genetic_code = 11)
observed <- translate_ggchord_cds(cds, genome, initiation = 'start')
expected <- read_fasta_lengths(
  'examples/format-inputs/protein-pair.faa', keep_sequence = TRUE
) |>
  filter(accver == 'BVI_023A_0009')
stopifnot(sub('[*]$', '', observed$models$protein_sequence) == expected$sequence)
letters <- observed$nucleotides |>
  filter(alignment_column <= 60)
protein <- observed$residues |>
  filter(alignment_column <= 20)
character_plot <- ggplot2::ggplot(
  letters, aes(alignment_column, accver, label = symbol)
) +
  geom_residue(aes(fill = symbol), size = 2.7, overflow = 'error') +
  ggplot2::geom_text(
    data = protein,
    aes(x = alignment_column * 3 - 1, y = 'translation', label = symbol),
    inherit.aes = FALSE, family = 'mono', size = 2.7
  ) +
  ggplot2::scale_fill_manual(
    values = c(A = '#B8D8BA', C = '#A6CBE3', G = '#F2D28C', T = '#EAA8A8'),
    guide = 'none'
  ) +
  ggplot2::scale_y_discrete(limits = c(coding$feature_id, 'translation')) +
  ggplot2::labs(x = 'Observed coding base (not genomic position)', y = NULL) +
  ggplot2::theme_minimal()
character_plot <- register_ggchord_alignment(character_plot, letters)
region <- letters |>
  summarise(accver = first(source_accver),
    start = min(source_position), end = max(source_position))
whole <- ggchord(genome |> select(accver, length), feature_data = coding) +
  geom_seq() +
  geom_gene() +
  coord_collinear()
result <- magnify_ggchord(
  whole, region, character_plot,
  center_size = c(.22, .7), detail_size = c(.5, .8)
)
# Full genomic codon membership, including negative-strand order:
protein |> select(alignment_column, symbol, codon, source_positions)
export_ggchord_fusion(result, width = 16, height = 5)$focus_regions
# ggplot2::ggsave('cds-detail.png', result, width = 16, height = 5)
