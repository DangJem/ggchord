# Run from the source checkout; optional dependencies: dplyr and Biostrings.
# Unchanged source records and provenance: examples/format-inputs/README.md.
library(ggchord)
library(dplyr)

genome <- read_fasta_lengths(
  'examples/format-inputs/emales-BVI_023A.fna', keep_sequence = TRUE
)
coding <- read_gff3('examples/format-inputs/emales-BVI_023A.gff3') |>
  filter(feature_id == 'BVI_023A_0009')
cds <- build_ggchord_cds(coding, genetic_code = 11)
translation <- translate_ggchord_cds(cds, genome, initiation = 'start')
letters <- translation$nucleotides |>
  filter(alignment_column <= 60)
region <- letters |>
  summarise(
    accver = first(source_accver),
    start = min(source_position), end = max(source_position)
  )
local_region <- coding |>
  select(accver, start, end)
focused <- focus_ggchord_data(genome, loci = local_region)
local <- ggchord(focused$seq_data) +
  geom_seq(show.legend = FALSE) +
  coord_collinear()
character <- ggplot2::ggplot(
  letters, aes(alignment_column, accver, label = symbol)
) +
  geom_residue(aes(fill = symbol), overflow = 'error', size = 2.7,
    show.legend = FALSE) +
  ggplot2::scale_fill_manual(
    values = c(A = '#B8D8BA', C = '#A6CBE3', G = '#F2D28C', T = '#EAA8A8')
  ) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    x = 'Observed coding base', y = NULL,
    title = '60 real negative-strand CDS bases'
  )
character <- register_ggchord_alignment(character, letters)
inner <- magnify_ggchord(
  local, region, character,
  center_box = c(.02, .2, .2, .8), detail_box = c(.27, .99, .15, .85)
)
whole <- ggchord(
  genome |> select(accver, length), feature_data = coding
) +
  geom_seq(show.legend = FALSE) +
  geom_gene() +
  coord_collinear()
outer_region <- local_region |>
  mutate(detail_panel = 'detail/center', purpose = 'Complete source CDS')
result <- magnify_ggchord(
  whole, outer_region, inner,
  center_box = c(.02, .22, .2, .8), detail_box = c(.27, .99, .1, .9)
)
report <- export_ggchord_fusion(result, width = 16, height = 6)
report$focus_regions
report$focus_region_data_by_composition
# ggplot2::ggsave('nested-cds.png', result, width = 16, height = 6, dpi = 120)
