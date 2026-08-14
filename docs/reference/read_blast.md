# Read one or more BLAST tabular output files into ribbon_data format

\`file\` reads a single BLAST tabular output file; \`files\` reads
multiple files at once and combines them. \`files\` accepts a character
vector of literal paths and/or wildcard patterns (e.g.
\`"examples/blastn/\*.o7"\`).

## Usage

``` r
read_blast(
  file = NULL,
  files = NULL,
  format = c("auto", "outfmt6", "outfmt7", "custom"),
  col_names = NULL,
  comment = "#",
  ...
)
```

## Arguments

- file:

  Optional path to a single BLAST tabular output file.

- files:

  Optional character vector of BLAST tabular output files (literal paths
  and/or wildcard patterns). All matched files are read and combined
  into one data.frame.

- format:

  Character. \`"auto"\` (default) detects the column layout from the
  number of columns; \`"outfmt6"\` / \`"outfmt7"\` require the standard
  12/17-column layouts; \`"custom"\` requires \`col_names\`.

- col_names:

  Optional character vector naming the columns in the file, used with
  \`format = "custom"\` or to override auto-detection.

- comment:

  Character comment character, default \`"#"\` (BLAST outfmt 7 header
  lines start with \`#\`).

- ...:

  Additional arguments passed to \[utils::read.table()\] (e.g.
  \`na.strings\`).

## Value

A data.frame with the required ribbon columns first (\`qaccver\`,
\`saccver\`, \`length\`, \`pident\`, \`qstart\`, \`qend\`, \`sstart\`,
\`send\`) followed by any preserved optional columns.

## Examples

``` r
library(ggchord)
blast_file <- tempfile(fileext = ".o7")
writeLines(c(
  "# BLASTN 2.13.0+",
  "seqA  seqB  98.5  1200  18  0  1  1200  1  1200  1e-180  2400"
), blast_file)
read_blast(blast_file)
#>   qaccver saccver length pident qstart qend sstart send mismatch gapopen evalue
#> 1    seqA    seqB   1200   98.5      1 1200      1 1200       18       0 1e-180
#>   bitscore
#> 1     2400
```
