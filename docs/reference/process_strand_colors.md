# Process gene color parameters in strand mode

Standardizes gene color parameters in strand mode (color by strand
direction) into a named vector with "+"/"-"

## Usage

``` r
process_strand_colors(gene_colors)
```

## Arguments

- gene_colors:

  Color vector (can be NULL, single value, vector of length 2, named
  vector with "+"/"-")

## Value

Named vector (names are "+"/"-"), standardized color values (default "+"
is red, "-" is blue)
