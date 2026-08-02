# Add arrowhead annotations at the tip of every sequence arc

plotly's scatter traces cannot draw line arrowheads, so the directional
arrows used by
[`geom_seq()`](https://dangjem.github.io/ggchord/reference/geom_seq.md)
are reproduced as plotly annotations.

## Usage

``` r
ggchord_plotly_arrows(pl, layout, colors = NULL)
```
