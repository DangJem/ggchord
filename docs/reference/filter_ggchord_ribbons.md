# Filter alignment ribbons before plotting

Keeps the ribbon rows that satisfy all requested criteria and optionally
sorts them. The returned data frame keeps every extra column and the
original column order; the original row numbers are attached as the
`source_rows` attribute.

## Usage

``` r
filter_ggchord_ribbons(
  ribbon_data,
  seq_ids = NULL,
  min_pident = NULL,
  max_pident = NULL,
  min_length = NULL,
  max_evalue = NULL,
  min_bitscore = NULL,
  min_query_coverage = NULL,
  min_subject_coverage = NULL,
  keep_pairs = NULL,
  drop_self_links = TRUE,
  sort_by = NULL
)
```

## Arguments

- ribbon_data:

  data.frame with at least `qaccver` and `saccver` columns (normally the
  full ribbon_data format).

- seq_ids:

  Optional character vector. Keep only rows where both the query and the
  subject are in `seq_ids`.

- min_pident, max_pident:

  Optional numeric. Lower/upper bounds on `pident`.

- min_length:

  Optional numeric. Lower bound on `length`.

- max_evalue:

  Optional numeric. Upper bound on `evalue` (requires an `evalue`
  column).

- min_bitscore:

  Optional numeric. Lower bound on `bitscore` (requires a `bitscore`
  column).

- min_query_coverage:

  Optional numeric (0-100). Lower bound on the query coverage. Uses the
  `qcovs` column when present; otherwise computed from `qlen` and the
  query interval (requires `qlen`).

- min_subject_coverage:

  Optional numeric (0-100). Lower bound on the subject coverage,
  computed from `slen` and the subject interval (requires an `slen`
  column).

- keep_pairs:

  Optional data.frame/list/matrix describing an undirected set of
  sequence pairs. A data.frame or matrix with the first two columns used
  as query/subject IDs, or a list of length-2 character vectors. A row
  is kept when its query/subject pair matches any pair in either
  direction.

- drop_self_links:

  Logical, default `TRUE`. Remove rows where `qaccver == saccver`.

- sort_by:

  Optional character vector of column names. Prefix a name with `"-"` or
  `"desc:"` for descending order (e.g. `c("pident", "-evalue")`).

## Value

A list with `data` (the filtered data frame, with `source_rows`
attribute) and `report` (n_input/n_kept/n_removed, removed_by_reason,
removed_rows and kept_rows).

## Examples

``` r
library(ggchord)
data(ribbon_data_example)
out <- filter_ggchord_ribbons(
  ribbon_data_example,
  min_pident = 95,
  drop_self_links = TRUE,
  sort_by = c("pident", "-length")
)
out$report
#> $n_input
#> [1] 31
#> 
#> $n_kept
#> [1] 3
#> 
#> $n_removed
#> [1] 28
#> 
#> $removed_by_reason
#>       reason  n
#> 1 min_pident 28
#> 
#> $removed_rows
#>  [1]  3  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28
#> [26] 29 30 31
#> 
#> $kept_rows
#> [1] 1 2 4
#> 
#> $sort_by
#> [1] "pident"  "-length"
#> 
```
