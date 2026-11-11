## Source tables and display copies {#source}

This chapter uses bundled sequence lengths, observed HSPs and source-strand
annotations. `feature_data` is a common annotation source. Generic feature layers
inherit the whole table; default gene layers require explicit gene/CDS classes.
An explicit layer table or data function takes priority. Use disjoint subsets
when combining the two kinds of geometry.

{{code:foundationdata}}

Prepare, clean, focus and segment helpers also accept `feature_data` and return
that member. Do not combine it with `gene_data`. Unknown strands remain unknown;
the legacy gene input retains its directional semantics and positional arguments.

## Named A–B–A routes {#instances}

The three displayed rows refer to two biological sequences. `seq_order` values
identify sources; names identify instances. Filter each layer's original HSPs
before routing. `between` accepts two names or two declaration-order indices.
Ambiguous routes fail. The two subsets below are disjoint and preserve every HSP.

{{code:foundationinstances}}

{{figure:foundationvertical}}

The same `foundation_plot` can use `coord_chord()` or `coord_collinear()`.
Vertical Collinear changes the global view while keeping source strands and
positions. Positive local normal points left. Explicit text styles retain
priority. `map_ggchord_position()` accepts an instance name; its inverse input
is final y in vertical views.

`export_ggchord_layout()` supplies `display_instances`, `instance_routing` and
`instance_omissions`. Explicit `source_hsp_id` detects repeated display across
layers; use `allow_repeated=TRUE` only when that repetition is intended. Without
an ID, identity stays layer-local; equal row values do not establish identity.

## Mirrored two-column relationships {#sankey}

Sankey displays each source once in each column. Query endpoints connect from
the left to subject endpoints on the right. Reciprocal input rows stay distinct;
same-entity links are excluded by default (`self_links=TRUE` opts in). Widths
come from both real intervals. This layout expresses relationships, not flow
conservation. It owns its left/right copies and does not accept `between` or
shared branches in this first version.

{{code:foundationsankey}}

{{figure:foundationsankey}}

## Count observations or display supplied values {#tracks}

`bar` counts observations at source positions; `col` displays supplied values.
Each accession/window computes its Stat and data Position before view projection
or display replication. `track=position_track(...)` selects the physical band;
`position` independently selects a native data adjustment such as stack or dodge.
The earlier `position=position_track(...)` band spelling remains accepted.

{{code:foundationcount}}

{{figure:foundationcount}}

The shortcut below counts the same observations using the same physical band.

```r
p_count_shortcut <- ggchord(foundation_seqs) +
  geom_seq() +
  geom_track_bar(
    aes(x = pos), data = foundation_events,
    track = position_track("query_starts")
  ) +
  coord_collinear()
```


| Geometry | Default Stat | Input meaning |
| --- | --- | --- |
| bar | count | x/pos observations, optional weight |
| col | identity | supplied y/track_value; interval or position |
| point, line, path, step, text | identity | source x and local y; text requires label |
| tile, rect, area, ribbon, segment, linerange, errorbar | identity | corresponding native ggplot2 extents |

Available statistics are identity/count/bin/summary. Bin needs explicit binwidth
or breaks. Native group and after_stat mappings are supported; numerical axes
retain the statistic's value units. Arbitrary third-party Geom/Stat plugins and
smooth/density/ecdf statistics are outside this first support matrix. Existing
track line/tile shortcuts retain their documented behaviour.

Source position uses x orientation; native text size currently uses mm and polygon outlines are full. Other orientations/size units/partial outlines fail explicitly. Native text styles and path arrows are forwarded.

Native tracks currently require one panel and a ggchord coordinate system. Facets, `coord_flip()`, raster, boxplot and violin are outside this first supported set and are rejected explicitly. Use `coord_collinear(seq_axis = "vertical")` for vertical projection. Statistical source membership is preserved before dodge or other data positions. Rectangles and tiles use their supplied bounds without an artificial zero baseline.

## Visual length and connection style {#projection}

`length_mode="equal"` changes angular allocation or linear span. Original lengths,
intervals and local axes remain unchanged. `sequence_display` exports both source
length and projection factor. A shared absolute scale bar is rejected for equal
multi-sequence displays. Default axes show 1 and the final source position;
intermediate ticks use available physical space. Closed circular seams combine
the two endpoints. Explicit breaks and labels take priority.

{{code:foundationequal}}

{{figure:foundationequal}}

Default ribbons keep Identity fills and full interval fronts.
Use `geom_link_ribbon(gradient = TRUE)` to interpolate endpoint sequence colours along the curve. `tip="subject"` opts into a decorative taper; its
shoulder retains the translated target span and its marker is not a single
matched base or a biological transfer direction. The source endpoints stay in layout exports.
Use `tip="none"` to keep full fronts. Fixed fill, mapped ribbon_fill or an explicit
ribbon fill scale overrides endpoint interpolation. Explicit pident uses a fixed
0–100 domain, boundaries 0/50/80/90/95/100 and five numerical interval proportions
50/30/10/5/5%. User scales and guides remain authoritative.

### Optional subject taper {#tip}

{{code:foundationarrow}}

{{figure:foundationarrow}}

Automatic categorical guides measure text to choose rows/columns, while the
render loop measures guide and panel space together. Explicit guide/theme
settings win. Circular outer labels share a measured queue at every label count;
`label_visibility` reports source row, visibility and available placement fields.
An unavailable collision count is not zero. Tiny canvases can remain crowded;
inspect [output diagnostics](output.html) at the intended export size.

Automatic Identity rulers reserve length from measured boundary text rather than generic legend swatch sizes. Set `legend.ribbon.key.height` or `legend.ribbon.key.width`, or supply an explicit guide, for a manual length. In circular layouts, parallel axis labels follow the circle glyph by glyph; horizontal and numeric orientations remain straight.
