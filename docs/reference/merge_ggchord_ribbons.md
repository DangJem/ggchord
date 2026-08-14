# Merge adjacent or overlapping alignment blocks of the same sequence pair

Merges alignment blocks that belong to the same (query, subject) pair,
are adjacent (gap \<= `max_gap`) or overlapping on both sequences, and
are compatible. The merged pident is weighted by alignment length.
Merging is deliberately conservative: blocks whose merged query and
subject spans would be inconsistent (unequal), whose pident differs by
more than `min_pident_difference`, or whose orientation differs (when
`require_same_orientation`) are left unmerged.

## Usage

``` r
merge_ggchord_ribbons(
  ribbon_data,
  max_gap = 0,
  min_pident_difference = 0,
  require_same_orientation = TRUE,
  group_by = c("qaccver", "saccver")
)
```

## Arguments

- ribbon_data:

  data.frame in ribbon_data format.

- max_gap:

  Numeric, default 0. Maximum gap (in bp) allowed between two blocks on
  both sequences for them to be merged.

- min_pident_difference:

  Numeric, default 0. When \> 0, two blocks are only merged when their
  pident difference is \<= this value.

- require_same_orientation:

  Logical, default `TRUE`. Only merge blocks whose query/subject
  direction is the same (both ascending or both descending, i.e.
  collinear or both inverted in the same way).

- group_by:

  Character vector, default `c("qaccver", "saccver")`. Columns used to
  identify the same sequence pair.

## Value

A list with `data` (merged data frame with `source_rows` attribute) and
`report` (data.frame with `output_row`, `from_rows` and `n_merged` for
every output row).

## Examples

``` r
# Two adjacent, collinear blocks of the same pair -> one merged block
rb <- data.frame(
  qaccver = c("A", "A"), saccver = c("B", "B"),
  length = c(100, 100), pident = c(95, 97),
  qstart = c(1, 101), qend = c(100, 200),
  sstart = c(501, 601), send = c(600, 700)
)
out <- merge_ggchord_ribbons(rb, max_gap = 0)
out$data
#>   qaccver saccver length pident qstart qend sstart send
#> 1       A       B    200     96      1  200    501  700
out$report
#>   output_row from_rows n_merged
#> 1          1       1,2        2
```
