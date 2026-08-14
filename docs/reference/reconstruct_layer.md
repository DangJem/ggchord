# Reconstruct a layer with the given data (and optional remapped mapping).

LayerInstance objects cannot be cloned with `ggproto(NULL, .)`, so the
layer is rebuilt through `layer()` with the same
geom/stat/mapping/params.

## Usage

``` r
reconstruct_layer(lyr, data, mapping = NULL)
```
