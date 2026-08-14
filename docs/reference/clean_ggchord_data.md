# Clean ggchord input data with explicit, report-driven policies

Applies a conservative, predictable set of cleaning policies to the
three ggchord tables and returns the cleaned copies plus a full report.
The original data frames are never modified. Nothing is dropped
silently: every change (including every dropped row) is recorded in
`report` with the original row number, the reason, the original value(s)
and the new value(s).

## Usage

``` r
clean_ggchord_data(
  seq_data,
  ribbon_data = NULL,
  gene_data = NULL,
  unknown_id = c("drop", "error", "keep"),
  out_of_range = c("clip", "drop", "error", "keep"),
  reversed_interval = c("sort", "drop", "error", "keep"),
  invalid_pident = c("clip", "drop", "error", "keep"),
  empty_annotation = c("keep", "drop", "replace"),
  replacement_annotation = "unannotated"
)
```

## Arguments

- seq_data:

  data.frame/tibble, required. Must contain `seq_id` and `length`; used
  as the coordinate reference for clipping.

- ribbon_data:

  data.frame/tibble, optional. Alignment results.

- gene_data:

  data.frame/tibble, optional. Gene annotation data.

- unknown_id:

  Character, default `"drop"`. Policy for rows whose sequence ID is
  missing, empty or unknown: `"drop"` removes them, `"error"` stops with
  a message, `"keep"` leaves them unchanged.

- out_of_range:

  Character, default `"clip"`. Policy for coordinates outside
  `[1, sequence length]`: `"clip"` clamps them, `"drop"` removes the
  row, `"error"` stops, `"keep"` leaves them unchanged. After clipping,
  intervals that become invalid (degenerate) are removed explicitly and
  reported.

- reversed_interval:

  Character, default `"sort"`. Policy for intervals with `start > end`:
  `"sort"` swaps the endpoints so the feature draws stably (the original
  direction is recorded in the report), `"drop"` removes the row,
  `"error"` stops, `"keep"` leaves them unchanged.

- invalid_pident:

  Character, default `"clip"`. Policy for `pident` outside `[0, 100]`:
  `"clip"` clamps them, `"drop"` removes the row, `"error"` stops,
  `"keep"` leaves them unchanged.

- empty_annotation:

  Character, default `"keep"`. Policy for missing/empty `anno` in gene
  data: `"replace"` fills them with `replacement_annotation`, `"drop"`
  removes the row, `"keep"` leaves them unchanged.

- replacement_annotation:

  Character, default `"unannotated"`. Replacement annotation used when
  `empty_annotation = "replace"`.

## Value

A list with four components: `seq_data`, `ribbon_data`, `gene_data`
(cleaned copies) and `report` (a data.frame with columns `table`, `row`
(original row number), `column`, `reason`, `original_value`, `new_value`
and `action`).

## Examples

``` r
data(seq_data_example)
data(ribbon_data_example)
data(gene_data_example)

# Introduce a few typical problems
bad_r <- transform(ribbon_data_example,
                   qstart = pmin(qstart, 1),
                   pident = pmin(pident, 150))
bad_g <- transform(gene_data_example,
                   anno = ifelse(seq_len(nrow(gene_data_example)) == 1,
                                 NA_character_, anno))

out <- clean_ggchord_data(seq_data_example, bad_r, bad_g)
head(out$report)
#>   table row column                             reason original_value new_value
#> 1  gene   1   anno missing or empty annotation (kept)           <NA>      <NA>
#>   action
#> 1   keep
# The cleaned tables are ready for ggchord()
# \donttest{
p <- ggchord(out$seq_data, out$ribbon_data, out$gene_data) +
  geom_seq() + geom_ribbon() + geom_gene()
#> Warning: ggchord(): input data has 1 validation warning(s) (e.g. "gene_data$anno is missing or empty (clean_ggchord_data(empty_annotation = 'replace') fills it)"). Run validate_ggchord_data(...) for details.
# }
```
