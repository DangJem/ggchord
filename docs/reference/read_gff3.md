# Read one or more GFF3 files into gene_data format

Parses a GFF3 file (9 tab-separated columns) and returns a data frame
that can be passed to
[`ggchord()`](https://dangjem.github.io/ggchord/reference/ggchord.md) as
`gene_data`. Only features whose `type` is in `feature_types` are kept.
The annotation label (`anno`) is extracted from the GFF3 attributes
column by trying the keys in `anno_from` in order (e.g. `product`, then
`Name`). The original useful columns (`type`, `source`, `score`,
`phase`, `attributes`) are preserved.

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
  patterns such as `"examples/gff3/*.gff3"`). All matched files are read
  and combined into one data.frame.

- feature_types:

  Character vector of feature types to keep, default `"CDS"`. Common
  choices: `"gene"`, `"tRNA"`, `"rRNA"`, `"repeat_region"`.

- anno_from:

  Character vector of GFF3 attribute keys, tried in order to fill the
  `anno` column, default `c("product", "Name", "gene", "ID")`. `anno` is
  `NA` when none of the keys is present.

- unstranded:

  Character, default `"plus"`. How to treat features whose strand is
  `"."` or `"?"`: `"plus"` maps them to `"+"`, `"drop"` removes them
  (ggchord requires `"+"` or `"-"`).

## Value

A data.frame with `seq_id`, `start`, `end`, `strand`, `anno` followed by
`type`, `source`, `score`, `phase` and `attributes`.

## Examples

``` r
# \donttest{
gff <- tempfile(fileext = ".gff3")
writeLines(c(
  "##gff-version 3",
  "seqA\tsource\tCDS\t101\t500\t.\t+\t0\tID=cds1;product=hypothetical protein",
  "seqA\tsource\tCDS\t600\t900\t.\t-\t0\tID=cds2;Name=integrase",
  "seqA\tsource\ttRNA\t1000\t1080\t.\t+\t0\tID=trna1"
), gff)
gd <- read_gff3(gff)
gd
#>   seq_id start end strand                 anno type source score phase
#> 1   seqA   101 500      + hypothetical protein  CDS source     .     0
#> 2   seqA   600 900      -            integrase  CDS source     .     0
#>                             attributes
#> 1 ID=cds1;product=hypothetical protein
#> 2               ID=cds2;Name=integrase
# }
```
