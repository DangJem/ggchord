# Deduplicate alignment ribbons

Removes fully duplicated, coordinate-near-duplicated or highly
overlapping alignment blocks within each (query, subject) pair and keeps
the best representative per `keep`.

## Usage

``` r
deduplicate_ggchord_ribbons(
  ribbon_data,
  tolerance = 0,
  by = c("exact", "coordinates", "overlap"),
  keep = c("best_pident", "longest", "first"),
  min_reciprocal_overlap = 0.9
)
```

## Arguments

- ribbon_data:

  data.frame in ribbon_data format.

- tolerance:

  Numeric, default 0. Maximum absolute difference (in bp) allowed on
  each of `qstart`/`qend`/`sstart`/`send` for two rows to count as
  coordinate duplicates (used with `by = "coordinates"`).

- by:

  Character, default `"exact"`. `"exact"` removes rows with identical
  pair and coordinates; `"coordinates"` additionally treats rows within
  `tolerance` bp as duplicates; `"overlap"` treats blocks whose
  reciprocal overlap is at least `min_reciprocal_overlap` on both
  sequences as duplicates.

- keep:

  Character, default `"best_pident"`. Which representative to keep:
  `"best_pident"` (highest pident), `"longest"` (longest length) or
  `"first"` (first occurrence in the input).

- min_reciprocal_overlap:

  Numeric (0-1), default 0.9. Reciprocal overlap threshold used with
  `by = "overlap"`.

## Value

A list with `data` (deduplicated data frame with `source_rows`
attribute) and `report` (n_input/n_kept/n_removed plus a data.frame of
removed rows with `row`, `duplicate_of` and `reason`).

## Examples

``` r
library(ggchord)
data(ribbon_data_example)
dup <- rbind(ribbon_data_example, ribbon_data_example[1, ])
out <- deduplicate_ggchord_ribbons(dup, by = "exact")
out$report
#> $n_input
#> [1] 32
#> 
#> $n_kept
#> [1] 31
#> 
#> $n_removed
#> [1] 1
#> 
#> $removed
#>   row duplicate_of          reason
#> 1  32            1 duplicate_exact
#> 
```
