# Read one or more GFF3 files into gene_data format

\`file\` reads a single GFF3 file; \`files\` reads multiple files at
once and combines them. \`files\` accepts a character vector of literal
paths and/or wildcard patterns (e.g. \`"examples/gff3/\*.gff3"\`).

## Usage

``` r
read_gff3(
  file = NULL,
  files = NULL,
  feature_types = "CDS",
  anno_from = c("product", "Name", "gene", "ID"),
  unstranded = c("plus", "drop")
)
```

## Arguments

- file:

  Optional path to a single GFF3 file.

- files:

  Optional character vector of GFF3 files (literal paths and/or wildcard
  patterns). All matched files are read and combined.

- feature_types:

  Character vector of feature types to keep, default \`"CDS"\`.

- anno_from:

  Character vector of GFF3 attribute keys, tried in order to fill the
  \`anno\` column.

- unstranded:

  Character, default \`"plus"\`.

## Value

A data.frame with \`seq_id\`, \`start\`, \`end\`, \`strand\`, \`anno\`
followed by \`type\`, \`source\`, \`score\`, \`phase\` and
\`attributes\`.

## Examples

``` r
library(ggchord)
gff <- tempfile(fileext = ".gff3")
writeLines(c(
  "##gff-version 3",
  "seqA  source  CDS  101  500  .  +  0  ID=cds1;product=hypothetical protein"
), gff)
read_gff3(gff)
#>   seq_id start end strand                 anno type source score phase
#> 1   seqA   101 500      + hypothetical protein  CDS source     .     0
#>                             attributes
#> 1 ID=cds1;product=hypothetical protein
```
