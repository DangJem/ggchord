<div class="task-intro"><p>Give colours a clear meaning: identify categories or continuous quantities, then choose mappings, scales and guides. Website accents do not dictate data colours.</p></div>

<span id="prepare"></span>

## Build styles from data roles {#roles}

Use four chord backbones and one observed alignment. Map the ribbon fill to the query ID and assign that category one named colour. Sequence colours retain their default scale.

{{code:styling}}
{{figure:styling}}

| Role | Mapped data | Useful scale family |
| --- | --- | --- |
| `seq_colour`, `seq_ring` | Sequence/category | Sequence colour or ring manual scales |
| `ribbon_fill`, `ribbon_alpha`, `ribbon_colour`, `ribbon_linetype` | Alignment metadata or statistics | Manual, gradient/steps/identity and alpha scales |
| `feature_fill`, `feature_shape`, `gene_fill`, `region_fill` | Annotation category or strand | Role-specific manual/preset scales |
| `link_colour`, `link_alpha`, `link_linewidth`, `link_linetype` | Point-link metadata | Continuous or manual link scales |
| `track_value`, `track_fill`, `track_colour` | Values/categories in one named track | Named value domain and colour/fill scales |

## Fixed styles and user scale priority {#priority}
Put a constant style outside `aes()`. Map variable styles inside it and supply an
appropriate scale. An explicitly fixed geom value takes precedence; omitting that
argument keeps mapping active. Your scale takes priority over automatic defaults.
Repeated same-kind layers keep independent data/mapping/parameters. A geometric
vector named by `accver` does not mean a scale also accepts per-sequence lists.

## Make legends readable {#guides}
`guide_ggchord_legend()` controls discrete role legends; `guide_ggchord_colourbar()`
controls continuous colour guides. The colour/color spelling aliases are retained,
with `colour` preferred in this documentation. Use `show.legend` on a layer,
`guides()` for guide choices and `theme(legend.position = ...)` for placement.
A wide guide may need multiple rows or a wider export, not smaller biological data.

## Select a theme for the scene {#themes}
This page uses the default theme. Add `theme()` to adjust text, legends or background; for example, `theme(text = element_text(size = 12))` sets the base text size. `theme_ggchord_publication()` provides a publication style and `theme_ggchord_plasmid()` suits circular maps. Themes control non-data appearance without changing coordinates or mapped colours.

## Target-size layout {#adaptive}

Add `adapt_ggchord()` to opt into target-size layout. The bundled example uses
recorded sequence lengths and HSPs; the reference canvas is in inches and the
automatic text floor is in mm. The report describes the requested 8 × 6 inch
PNG; use the same size with standard `ggsave()` to reproduce that arrangement.

{{code:courseblockw}}

`fixed` holds default theme and guide furniture at reference physical sizes.
`scale` adjusts theme-driven defaults within the chosen multiplier range.
`adaptive` adds a bounded joint search for native guide grids, sides and
directions, then remeasures the result. It reuses existing label layout;
`feasible` is a guide-space result, not a promise that every label is collision-free.
No policy preserves the existing plotting behaviour.

Explicit geom styles, theme sizes, guide positions/directions/grids/key sizes,
scales and `show.legend` remain authoritative. Continuous/stepped bars remain
intact; categories and source records are retained. When constraints exceed the
canvas budget, diagnostics report `insufficient_space` and rendering warns.
Increase the export size or revise the explicit guide settings. The bounded
search does not guarantee a global optimum. Layout-only exports record the
policy but require a target device for final guide measurements.

Apply the policy to each desired Fusion panel before composing it. Each panel
uses its allocated size; outer panel allocations and the cross-panel guide are
separate. Facets, arbitrary custom grobs, blanket scaling of geometric constants
and font embedding certification are outside this policy.

## Result checks {#check-result}
Links from the same query share a colour; continuous identity needs a continuous scale. Put fixed colours outside aes() and make each guide explain its mapping.

<div class="next-step">Continue: <a href="annotations.html">Features, stacking and genomic decoration →</a></div>

See [scale families](scale-families.html) for runnable manual/continuous/steps/identity, line type, opacity and theme examples.

## Detailed appearance with standard theme syntax {#native-theme}

`ggchord()` already supplies `theme_ggchord()`. Choose one of the four presets only when its appearance suits the figure: Minimal reduces decoration; Publication uses compact typography; Dark changes non-data decoration for a dark background; Plasmid emphasises circular annotations. Palettes, scientific relationships and sequence directions remain independent. Use explicit data scales or fixed geom colours when a background requires different data contrast.

| Change | Component |
| --- | --- |
| Background, font, non-data labels/leader appearance | `theme()` |
| Sequence colours, feature fill, numerical Track range/breaks/units | `aes()` and the corresponding scale |
| Backbone form, Ribbon arrow/gradient mode, fixed data styles | Geom |
| Track band, feature stacking and offsets | Position |
| Sequence order/direction and projection | Coord / explicit data workflow |
| Guide rows, columns, order, direction and individual position | Guide; global placement in theme |

The registered `ggchord.axis.*` physical settings accept scalar absolute units, just like the existing short arguments in `theme_ggchord()`. Prefer `theme(ggchord.axis.gap = unit(1, "mm"))` for later adjustments. Use one spelling in a theme constructor; supplying both spellings for the same setting is an error. A later `theme()` merges fields in the familiar ggplot2 order. `ggchord.seq.label`, `ggchord.gene.label`, `ggchord.feature.label` and their `.segment` line elements style their actual annotation layers, while explicit geom styles/mappings remain authoritative. Ordinary Cartesian `axis.text.x` is not a substitute for the automatic genomic sequence axis; Dotplot has independent native x/y axes.

### Shared Track appearance and local overrides {#track-theme}

The following real GC fractions come from the unchanged bundled pUC19 sequence. It needs optional dplyr. The Track's data scale owns the fraction domain and unit; its local native Guide owns axis appearance. The ordinary theme owns shared decorations. No coverage or read support is invented.

{{code:courseblockx}}

`ggchord.track.axis.line/text/title` control shared numerical-axis styles; `.text.offset` and `.title.offset` are nonnegative absolute units. Resolved font/family/angle metrics are used when reserving plot space. The GC Guide's explicit text style wins over the later shared 9pt default. A local `element_blank()` hides just that decoration. `guide_axis()` supports local `axis.line/text/title`, `angle` and `check.overlap`; other options fail explicitly instead of being silently ignored.

Layers with one Track ID share a band, value domain and one numerical reference; distinct IDs already have independent value/fill/colour scales and colour guides. Add `scale_track_*()` for each ID rather than a theme element named after every track. Colour-guide appearance can use the native `guide_legend(theme = ...)` or `guide_colourbar(theme = ...)` supplied to the corresponding Track colour scale. Global legend fonts/background/placement remain theme concerns. Sequence order, HSP selection, statistics and units never change when switching themes. Automatic label/Guide layout is a separate policy and remeasures at the final export size.

<div class="task-intro"><p>Put common adjustments in the appropriate ggplot2 component: mappings, geometry, positions, coordinates and themes have distinct roles.</p></div>

<span id="mapping-prepare"></span>

## Assign each responsibility {#grammar}
[Plotting concepts](modes.html#prepare) explains layers, mappings, scales and coordinates. The example below applies that distinction to alignment identity.

## Map link identity {#mapping}


{{code:adjustments}}

{{figure:mappingexample}}

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

## Result checks {#check-result-supplement}
Use a coordinate for global rotation, named values for per-sequence geometry, and mappings plus scales for data-dependent colour.
