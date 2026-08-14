# Validate ggchord input data before plotting

Performs structured validation of `seq_data`, `ribbon_data` and
`gene_data` so that problems can be found, understood and fixed before
plotting. The result is a `"ggchord_validation"` object with a `valid`
flag, `errors` (severe problems that make the plot misleading),
`warnings` (drawable but noteworthy issues), per-category counts, a data
summary, the original row numbers of every problem, and a list of
automatically fixable issues.

## Usage

``` r
validate_ggchord_data(
  seq_data,
  ribbon_data = NULL,
  gene_data = NULL,
  strict = FALSE,
  check_coordinates = TRUE,
  check_duplicates = TRUE,
  check_self_links = TRUE
)
```

## Arguments

- seq_data:

  data.frame/tibble, required. Basic sequence information (columns
  `seq_id`, `length`).

- ribbon_data:

  data.frame/tibble, optional. Alignment results (columns `qaccver`,
  `saccver`, `length`, `pident`, `qstart`, `qend`, `sstart`, `send`).

- gene_data:

  data.frame/tibble, optional. Gene annotation data (columns `seq_id`,
  `start`, `end`, `strand`, `anno`).

- strict:

  Logical. When `TRUE`, stop with an error as soon as any severe problem
  is found. When `FALSE` (default), return the full diagnostic report
  without stopping.

- check_coordinates:

  Logical, default `TRUE`. Whether to check that ribbon/gene coordinates
  stay inside `[1, sequence length]`.

- check_duplicates:

  Logical, default `TRUE`. Whether to look for fully duplicated,
  near-duplicated and highly overlapping records.

- check_self_links:

  Logical, default `TRUE`. Whether to flag alignment rows where
  `qaccver == saccver`.

## Value

A `"ggchord_validation"` object (a list) with at least:

- `valid`:

  Logical: `TRUE` when there are no severe errors.

- `errors`:

  data.frame of severe issues (table, category, row, column, message).

- `warnings`:

  data.frame of non-severe issues (same columns).

- `summary`:

  Per-category counts (table, category, severity, n).

- `data_summary`:

  Counts of sequences/ribbons/genes, unknown IDs, out-of-range rows,
  etc.

- `invalid_rows`:

  Named list of original row numbers per problem category.

- `cleanable`:

  data.frame of fixable issues with suggested actions.

[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) methods are provided.

## Examples

``` r
data(seq_data_example)
data(ribbon_data_example)
data(gene_data_example)

res <- validate_ggchord_data(seq_data_example, ribbon_data_example,
                             gene_data_example)
res$valid
#> [1] TRUE
print(res)
#> ggchord data validation
#> ========================
#> Result: VALID (0 error(s), 0 warning(s))
#> Sequences: 4 | Ribbons: 31 | Genes: 20
summary(res)
#> ggchord validation summary
#> Valid: TRUE
#> Sequences: 4 | Ribbons: 31 | Genes: 20
#> No issues found.

# Introduce a problem: unknown sequence ID in the ribbons
bad <- transform(ribbon_data_example, saccver = "NOT_A_SEQUENCE")
v <- validate_ggchord_data(seq_data_example, bad)
v$invalid_rows$ribbon_unknown_id
#> NULL
```
