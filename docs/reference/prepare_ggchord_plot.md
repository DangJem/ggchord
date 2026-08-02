# Fully prepare a ggchord plot and return it (compute layout, rename ribbon mappings, attach scales, set coordinates). The layout is cached on the plot (and on the shared reference environment) during preparation. Used by the lazy layer data path so that plotly::ggplotly() sees the same state as a normal build.

Fully prepare a ggchord plot and return it (compute layout, rename
ribbon mappings, attach scales, set coordinates). The layout is cached
on the plot (and on the shared reference environment) during
preparation. Used by the lazy layer data path so that plotly::ggplotly()
sees the same state as a normal build.

## Usage

``` r
prepare_ggchord_plot(plot)
```
