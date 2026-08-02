# Custom gene arrow legend drawing function

Generates gene arrow-shaped legend symbols (polygons) for ggplot2
legends

## Usage

``` r
draw_key_gene_arrow(data, params, size)
```

## Arguments

- data:

  Legend data (contains aesthetic mapping parameters like fill, colour,
  size)

- params:

  Legend parameters (automatically passed by ggplot2)

- size:

  Legend symbol size

## Value

grid::polygonGrob object, gene arrow-shaped legend symbol
