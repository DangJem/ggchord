# Build a concrete sequence-group label layer from computed layout data

Group labels are appended at build time (not when the user calls
[`geom_seq()`](https://dangjem.github.io/ggchord/reference/geom_seq.md)),
which keeps
[`geom_seq()`](https://dangjem.github.io/ggchord/reference/geom_seq.md)
backward compatible: it always returns a single sequence-arc layer.

## Usage

``` r
ggchord_group_label_layer(group_labels)
```
