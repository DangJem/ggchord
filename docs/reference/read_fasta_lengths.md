# Read one or more FASTA files into seq_data format

\`file\` reads a single FASTA file; \`files\` reads multiple files at
once and combines them. \`files\` accepts a character vector of literal
paths and/or wildcard patterns (e.g. \`"examples/fasta/\*.fna"\`).

## Usage

``` r
read_fasta_lengths(file = NULL, files = NULL, header_delim = NULL)
```

## Arguments

- file:

  Optional path to a single FASTA file.

- files:

  Optional character vector of FASTA files (literal paths and/or
  wildcard patterns). All matched files are read and combined.

- header_delim:

  Optional character. When given, each header is split at every
  occurrence of this delimiter and only the first piece is kept.

## Value

A data.frame with columns \`seq_id\` and \`length\`.

## Examples

``` r
library(ggchord)
fasta <- tempfile(fileext = ".fna")
writeLines(c(">seqA some description", "ACGTACGTACGTACGT", "ACGTACGT",
             ">seqB", "TTTTGGGG"), fasta)
read_fasta_lengths(fasta)
#>   seq_id length
#> 1   seqA     24
#> 2   seqB      8
```
