# Convert a ggchord plot to a plotly object

ggchord provides an S3 method for
[`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html) so
that chord diagrams built with
[`ggchord()`](https://dangjem.github.io/ggchord/reference/ggchord.md)
(including plots that combine the ribbon and the gene layers) can be
converted to interactive plotly charts. The conversion first renders the
computed geometry with standard ggplot2 layers whose colors are
pre-mapped, then delegates to
[`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html).

## Usage

``` r
# S3 method for class 'ggchord'
ggplotly(p, ...)
```

## Arguments

- p:

  A ggchord plot object created with
  [`ggchord()`](https://dangjem.github.io/ggchord/reference/ggchord.md).

- ...:

  Additional arguments passed to
  [`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html).

## Value

A plotly object.
