# Read one or more BLAST tabular output files into ribbon_data format

Parses BLAST `-outfmt 6` or `-outfmt 7` tabular output into a data frame
that can be passed directly to
[`ggchord()`](https://dangjem.github.io/ggchord/reference/ggchord.md) as
`ribbon_data`. The standard 12-column layout
(`qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore`)
and the 17-column extension used by the README workflow
(`...+ qcovs qlen slen sstrand stitle`) are recognised automatically.
Optional columns such as `evalue`, `bitscore`, `qcovs`, `qlen`, `slen`,
`sstrand` and `stitle` are preserved for filtering and hover
diagnostics.

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
  and/or wildcard patterns such as `"examples/blastn/*.o7"`). All
  matched files are read and combined into one data.frame.

- format:

  Character. `"auto"` (default) detects the column layout from the
  number of columns; `"outfmt6"` / `"outfmt7"` require the standard
  12/17-column layouts; `"custom"` requires `col_names`.

- col_names:

  Optional character vector naming the columns in the file, used with
  `format = "custom"` or to override auto-detection.

- comment:

  Character comment character, default `"#"` (BLAST outfmt 7 header
  lines start with `#`).

- ...:

  Additional arguments passed to
  [`read.table`](https://rdrr.io/r/utils/read.table.html) (e.g.
  `na.strings`).

## Value

A data.frame with the required ribbon columns first (`qaccver`,
`saccver`, `length`, `pident`, `qstart`, `qend`, `sstart`, `send`)
followed by any preserved optional columns.

## Examples

``` r
# \donttest{
blast_file <- tempfile(fileext = ".o7")
writeLines(c(
  "# BLASTN 2.13.0+",
  "# Query: seqA",
  "seqA\tseqB\t98.5\t1200\t18\t0\t1\t1200\t1\t1200\t1e-180\t2400\t100\t5000\t4800\tplus\tseqB",
  "seqA\tseqC\t95.0\t800\t40\t2\t1300\t2100\t50\t850\t1e-100\t1500\t80\t5000\t6000\tminus\tseqC"
), blast_file)
rb <- read_blast(blast_file)
head(rb)
#>   qaccver saccver length pident qstart qend sstart send mismatch gapopen evalue
#> 1    seqA    seqB   1200   98.5      1 1200      1 1200       18       0 1e-180
#> 2    seqA    seqC    800   95.0   1300 2100     50  850       40       2 1e-100
#>   bitscore qcovs qlen slen sstrand stitle
#> 1     2400   100 5000 4800    plus   seqB
#> 2     1500    80 5000 6000   minus   seqC
# }
```
