# Process axis label orientation parameters

Standardizes axis label orientation parameters in various formats
(character/numeric/vector) into a named vector (mapped by sequence ID)

## Usage

``` r
process_axis_orientation(param, seqs)
```

## Arguments

- param:

  Character ("horizontal", "parallel" or "perpendicular"), numeric
  (angle), vector (length matches number of sequences), or named vector,
  label orientation parameter

- seqs:

  Character vector, list of sequence IDs

## Value

Named vector (names are seq_id), values are "horizontal", "parallel",
"perpendicular" or numeric angles
