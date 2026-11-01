## Assign each responsibility {#grammar}
| Task | Interface |
| --- | --- |
| Map a column to an aesthetic | `aes()` |
| Define mapping and guide | `scale_*()`, `guides()` |
| Geometry and fixed style | `geom_*()` |
| Local placement | `position_*()` |
| Global coordinate system | `coord_*()` |
| Non-data appearance | `theme_*()`, `theme()` |
| Titles and captions | `labs()` |

Fixed styles may be supplied directly to geoms. Data-driven variation needs a
mapping and scale. An explicit fixed value overrides its mapping; omitting the
fixed argument preserves the mapping. User scales override inferred defaults.
`shape` is the fixed feature geometry parameter; `feature_shape` remains the
mapped role. Explicit `shape = NULL` is invalid, not a request to use a mapping.

## Map link identity {#mapping}
After the [shared setup](index.html#quick-start):

{{code:adjustments}}

## Place features {#position}
Use `position_strand(c("+" = .1, "-" = -.2))` for signed strand offsets.
`position_feature_stack()` handles overlapping features; `position_track()`
allocates a named quantitative lane. Geometric sequence parameters accept one
value, a matching unnamed vector, or a vector named by `accver`; named vectors
are recommended. This syntax does not extend to scales, guides or themes.

Each layer owns its own `data`, `mapping` and parameters. Repeated feature or
link layers do not replace one another. Use ordinary ggplot2 `show.legend` and
`inherit.aes`; role input mappings precede geometry, while staged visual mappings
follow ggplot2 evaluation. See the [API contract](reference/ggchord_api_contract.html).
