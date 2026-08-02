# Process gene color parameters in manual mode

Standardizes gene color parameters in manual mode (color by gene
annotation) into a vector named by gene annotation

## Usage

``` r
process_manual_colors(gene_colors, unique_anno, gene_order)
```

## Arguments

- gene_colors:

  Color vector (can be NULL, single value, vector, named vector with
  gene annotations)

- unique_anno:

  Character vector, unique gene annotation names

- gene_order:

  Character vector, display order of genes in the legend, default NULL
  (order of appearance)

## Value

Named vector (names are gene annotations), standardized color values
(default uses the built-in Set1 palette)
