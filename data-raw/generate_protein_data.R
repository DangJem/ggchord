# Rebuild bundled protein tables from complete sequences and observed BLASTp rows.
# Run prepare-bundled-proteins.py first to repeat sequence selection and BLASTp.
pkgload::load_all(quiet = TRUE)
library(dplyr)
fasta_file <- "examples/protein/four-proteins.faa"
protein_seq_data_example <- read_fasta_lengths(fasta_file, sequence_unit = "aa") |>
  select(accver, length, sequence_unit)
fasta <- readLines(fasta_file)
records <- split(fasta, cumsum(startsWith(fasta, ">")))
protein_seq_data_example$sequence <- vapply(records,
  function(record) paste(record[-1], collapse = ""), character(1))
stopifnot(all(nchar(protein_seq_data_example$sequence) == protein_seq_data_example$length))

protein_ribbon_data_example <- read_blast("examples/protein/all-vs-all-blastp.tsv") |>
  mutate(qrank = match(qaccver, protein_seq_data_example$accver),
         srank = match(saccver, protein_seq_data_example$accver)) |>
  filter(qrank < srank) |>
  select(-qrank, -srank) |>
  mutate(sequence_unit = "aa")
stopifnot(nrow(protein_seq_data_example) == 4L,
          nrow(protein_ribbon_data_example) == choose(4, 2))
stopifnot(validate_ggchord_data(protein_seq_data_example,
                              protein_ribbon_data_example)$valid)
save(protein_seq_data_example, file = "data/protein_seq_data_example.rda",
     compress = "xz", version = 2)
save(protein_ribbon_data_example, file = "data/protein_ribbon_data_example.rda",
     compress = "xz", version = 2)
