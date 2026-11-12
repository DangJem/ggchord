# Optional dependency: dplyr. Source: unchanged bundled pUC19 FASTA.
library(ggchord)
library(dplyr)
genome <- read_fasta_lengths('examples/plasmid/pUC19.fna', keep_sequence = TRUE)
source_length <- genome$length[[1L]]
region <- tibble(accver = genome$accver, start = source_length - 46,
  end = 35, wrap = TRUE)
# Extract actual source characters, keeping their order across the origin.
positions <- c(seq.int(region$start, source_length), seq_len(region$end))
fragment <- paste0(substring(genome$sequence, positions, positions), collapse = '')
source <- tibble(accver = 'origin-window', source_accver = genome$accver,
  start = region$start, strand = '+', sequence = genome$sequence)
cells <- as_ggchord_alignment(c('origin-window' = fragment), 'DNA',
  source_data = source, circular = TRUE)
characters <- ggplot2::ggplot(cells, aes(alignment_column, accver, label = symbol)) +
  geom_residue(aes(fill = symbol), show.legend = FALSE) +
  ggplot2::theme_void() +
  ggplot2::labs(title = 'pUC19: source origin within an 82-base window')
characters <- register_ggchord_alignment(characters, cells)
whole <- ggchord(genome) + geom_seq(show.legend = FALSE) + coord_circular()
origin_detail <- magnify_ggchord(whole, region, characters,
  center_box = c(.01, .28, .05, .95), detail_box = c(.32, .99, .1, .9))
origin_detail
