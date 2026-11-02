# Run from a directory containing GFF3 files. Source records stay unchanged.
library(dplyr)
library(readr)
library(stringr)

gff_columns <- c(
  "accver", "source", "type", "start", "end", "score", "strand", "phase", "attributes"
)
gff_files <- list.files(".", pattern = "\\.gff3$", full.names = TRUE)
if (!length(gff_files)) stop("No GFF3 files found in the current directory.")

gff_rows <- lapply(gff_files, function(path) {
  read_tsv(path, col_names = gff_columns, comment = "#", show_col_types = FALSE)
}) |>
  bind_rows()

gene_track <- gff_rows |>
  filter(type == "CDS") |>
  mutate(anno = str_extract(attributes, "(?<=product=)[^;]+")) |>
  select(accver, start, end, strand, anno)

write_tsv(gene_track, "gene_track.tsv")
