# Read one or more FASTA files into seq_data format

Reads a FASTA file and returns a `seq_id`/`length` data frame suitable
for
[`ggchord()`](https://dangjem.github.io/ggchord/reference/ggchord.md).
The sequence ID is the first whitespace-delimited token of each header
by default; pass `header_delim` (e.g. `"|"`) to split NCBI-style headers
further and keep only the first field.

## Usage

``` r
read_fasta_lengths(file = NULL, files = NULL, header_delim = NULL)
```

## Arguments

- file:

  Optional path to a single FASTA file.

- files:

  Optional character vector of FASTA files (literal paths and/or
  wildcard patterns such as `"examples/fasta/*.fna"`). All matched files
  are read and combined into one data.frame.

- header_delim:

  Optional character. When given, each header is split at every
  occurrence of this delimiter and only the first piece is kept (e.g.
  `"|"` for NCBI headers such as `>NC_000001.1|cds|...`).

## Value

A data.frame with columns `seq_id` and `length`. A warning is emitted
when the file contains duplicate sequence IDs.

## Examples

``` r
fasta <- tempfile(fileext = ".fna")
writeLines(c(
  ">seqA some description",
  "ACGTACGTACGTACGT",
  "ACGTACGT",
  ">seqB",
  "TTTTGGGG"
), fasta)
read_fasta_lengths(fasta)
#>   seq_id length
#> 1   seqA     24
#> 2   seqB      8
```
