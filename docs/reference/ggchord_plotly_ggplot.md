# Build a plotly-ready standard ggplot2 plot from a ggchord plot

The returned plot uses standard geoms with pre-mapped colors and
identity scales, so plotly::ggplotly() can convert every element
(sequence arcs, ribbons, gene arrows, axis) without scale conflicts.
User-supplied color scales (added with `+`) are honored.

## Usage

``` r
ggchord_plotly_ggplot(p)
```
