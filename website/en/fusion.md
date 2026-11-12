## Connect real source intervals {#intervals}

Fusion combines independent genomic plots and connects their real sequence
intervals. The center can be Circular, Chord or Collinear; peripherals are horizontal or
vertical Collinear plots. Each plot remains usable on its own. The composition
rebuilds all panels at their actual physical sizes before drawing links.

Download the [original EMALE PAF and provenance](formats.html). This lesson
requires the optional `dplyr` package. The PAF contains 11 observed HSPs; the
explicit display filter below selects two with length ≥900 and identity ≥85%.
The count track bins all original query-start observations in explicit 2,000 bp bins, not coverage. The last display bin is explicitly clipped to the true sequence boundary; source events and counts are retained.

{{code:fusiondata}}

Keep reverse endpoints as normalized by `read_paf()`. A constructor warning
about descending endpoints does not mean they should be sorted. Prepare two
independent plots; explicit sequence colours keep their accession meanings
consistent. Local legends are hidden explicitly because their sequence colours
are already visible and only the cross-panel Identity legend is needed here.

{{code:fusionpanels}}

{{code:fusionlinear}}

The fixed `qpanel` and `spanel` select the exact query and subject panels.
The ribbon is anchored to the matching sequence in the multi-sequence panel,
not its outer rectangle or center. Query/subject roles do not indicate a
biological flow of genetic material.

{{code:fusionview}}

{{figure:fusionview}}

## Register panels and instances {#panels}

Use one panel per side: `right`, `top`, `left`, `bottom`. IDs must be unique;
`center` is reserved. `size=c(width,height)` and `center_size` are fractions of
the composition area after the cross-panel legend has been allocated. Explicit
allocations that overlap fail rather than silently relocating a panel.
A side specifies placement; it does not rotate a plot's text. Set the plot's
own `coord_collinear(seq_axis="vertical")` for a vertical sequence axis.

The four-direction example reuses source PAF rows 3, 2, 9 and 1, explicitly
connecting them north, west, south and east. These repeated views are separate
display instances, not extra observations. Different panels need not share a
physical bp/mm scale.

{{code:fusionfour}}

{{figure:fusionfour}}

Selectors can also be columns `qpanel`, `spanel`, `q_instance`, `s_instance`.
Do not provide both a fixed panel argument and its column. If an accession
occurs in several panels or repeated display instances, omission is allowed
only when exactly one fully visible interval matches. Ambiguity is an error;
name the intended panel and instance. Original source coordinates also work
for explicitly renamed focus windows via their retained source origins.
Repeated `source_hsp_id` display requires `allow_repeated=TRUE` in every
participating layer. Without IDs, source rows remain layer-local; identical
values are not treated as a shared biological identity.

## Styles and physical layout {#styles}

Cross-panel layers accept `aes(fill=..., colour=..., alpha=..., linewidth=...,
linetype=...)` with ordinary ggplot2 scales. Panel role-specific scales remain
local. Explicit fixed styles override their mappings. Default Identity is
0–100 with 0/50/80/90/95/100 proportional breaks. Endpoint gradients and
subject tips are opt-in, using the same renderer as single-view ribbons.
A tip is a decorative direction marker; full source fronts remain exported.

{{code:fusionstyle}}

Default `layout="fixed"` allocates panel rectangles and rebuilds their labels
and guides at the output size. `layout="adaptive"` additionally fits the
central rectangle to the physical aspect ratio. Explicit font sizes remain
fixed. This is not a general font-scaling or global collision solver. Each
panel retains its own guides; different scale meanings are not automatically
merged. Use `show.legend=FALSE` or local `theme()` deliberately where appropriate.
Composition `theme()`/`guides()` affect the cross-panel legend only. Required
guide space is measured; a canvas too small for it fails with a diagnostic.

## Export and inspect anchors {#export}

{{code:fusionexport}}

`panels` records actual panel/container bounds in mm. `connections` retains
original relationship columns, source rows and resolved routes; unknown
directions are explicitly omitted and reported. `fronts` contains local
coordinates, source interval boundaries, global mm and local tangent/normal
vectors. `paths` adds the rendered connection geometry; `source_layouts`,
`source_data` and `relationship_data` preserve their separate provenance.
Diagnostics include local label visibility, fit results and guide dimensions;
unmeasured global collisions are not reported as zero.

Intervals are one-based inclusive, with full boundaries `[min(start,end)-1,
max(start,end)]` ordered by the original endpoint direction. A one-position
interval retains nonzero width. This new composition rule does not rewrite
legacy single-view ribbon geometry. Source positions, interval boundaries and
mm coordinates are different quantities.

`export_ggchord_fusion()` measures on a temporary device and restores the
previous device. Saving rebuilds on the actual output device, so font metrics
may differ slightly while source coordinates remain unchanged. Use ordinary
`ggsave()` or the shared preview; do not use `ggplot_build()` on a composition.

{{code:fusionsave}}

Fully visible intervals are required. Prepare partial windows or splits in
the source workflow explicitly. Cross-origin wrapping, direct bp↔aa links,
peripheral Sankey/Radial/Dotplot, facets, external ggplot coordinate guessing,
free-angle/eight-direction layouts and automatic global routing are outside
this interface. Links may cross or pass through panel interiors; foreground
annotations are retained, but arbitrary layouts are not guaranteed collision-free.


## A source interval and its local detail {#magnify}

`magnify_ggchord()` composes a Circular/Chord/Collinear source and an explicit Collinear detail plot. A source interval is biological data, not a rectangle in image pixels. The same region must be fully visible in both panels. Both rebuild at the final device size; no screenshot is stretched.

This example selects pETDuet-1 source positions 101–1900 (1800 bp). Its origin-crossing annotation at 5404→2 has no overlap with this interval; exclude it explicitly before linear focusing rather than taking its min/max as a long feature. The detail retains 10 annotations and labels the explicitly selected CDS/replication-origin subset. Reversing the detail's display orientation changes neither source coordinates nor biological strands.

```r
plasmid <- plasmid_example_pETDuet_1
annotations <- find_common_features(plasmid)
selected_region <- data.frame(accver = plasmid$accver, start = 101, end = 1900)
local_source <- focus_ggchord_data(
  plasmid,
  feature_data = annotations |> filter(start <= end),
  loci = selected_region
)
whole_plasmid <- ggchord(plasmid, feature_data = annotations) +
  geom_seq(show.legend = FALSE) +
  geom_feature(position = position_plasmid(), show.legend = FALSE) +
  geom_seq_center_label() +
  coord_circular(rotation = -30)
local_labels <- local_source$feature_data |>
  filter(type %in% c("CDS", "rep_origin"), anno != "ATG")
local_detail <- ggchord(
  local_source$seq_data, feature_data = local_source$feature_data
) +
  geom_seq(seq_orientation = -1, show.legend = FALSE) +
  geom_feature(show.legend = FALSE) +
  geom_feature_label(
    data = local_labels, label_orientation = "horizontal",
    label_overlap = "allow", colour = "#202020", size = 3
  ) +
  coord_collinear(seq_axis = "vertical") +
  scale_seq_position_continuous(coordinates = "source")
focused_view <- magnify_ggchord(whole_plasmid, selected_region, local_detail)
focused_view
focus_geometry <- export_ggchord_fusion(focused_view, width = 12, height = 8)
focus_geometry$focus_regions
focus_geometry$leaders
```

Use ordinary `ggsave()` at a chosen physical size, for example 12×8 inches; `view_ggchord()` uses the same renderer. `center_size` and `detail_size` are container fractions, while `highlight_width` is physical mm. Changing the detail side or orientation can change leader crossings; global collision avoidance is not claimed. `connector = "none"` keeps only highlights.

The dashed lines are decorative focus leaders. They are exported separately from biological `connections` and `relationship_data`; adding actual `fusion_ribbon()` layers still requires real relationship data. `focus_region_data` preserves the original selection and arbitrary metadata, while `focus_regions` reports resolved instances, source origins and the sampled physical path-length ratio. The ratio can be below one with a small detail allocation; it is not a change to biological length.

The detail plot owns its data, scales, layers and statistics. For count/bin/coverage tracks, prepare and state the intended local or global statistical population yourself; the wrapper records `explicit_detail_plot` and never silently chooses a scope. It does not copy arbitrary source layers, infer alignments or accept ordinary Cartesian character panels. Use one ordered region, explicitly select `instance`/`detail_instance` when an accession is repeated, and split origin wraps before this single-window workflow. Multiple/nested inset boxes and other central Coords remain outside this interface.

## Paired source windows with real HSPs {#paired-detail}

Each row of `region` selects one original source interval; all rows share one explicit detail panel. The example selects the two endpoints of a real HSP from the four-sequence dataset. Gene direction and the local HSP come from the focused input tables. `region_row` links each derived region, front and decorative leader to its input row; arbitrary original metadata stays in `focus_region_data`.

```r
library(dplyr)
seqs <- seq_data_example |> select(accver, length)
features <- gene_data_example |> mutate(type = 'gene')
hsp <- ribbon_data_example |> slice(1)
regions <- bind_rows(
  transmute(hsp, accver = qaccver, start = pmin(qstart, qend), end = pmax(qstart, qend)),
  transmute(hsp, accver = saccver, start = pmin(sstart, send), end = pmax(sstart, send))
)
local <- focus_ggchord_data(seqs, feature_data = features,
  ribbon_data = hsp, loci = regions)
source <- ggchord(seqs, feature_data = features, ribbon_data = hsp) +
  geom_seq(show.legend = FALSE) + geom_gene(show.legend = FALSE) +
  geom_link_ribbon(show.legend = FALSE)
detail <- ggchord(local$seq_data, feature_data = local$feature_data,
  ribbon_data = local$ribbon_data) + geom_seq(show.legend = FALSE) +
  geom_gene(show.legend = FALSE) + geom_link_ribbon(show.legend = FALSE) + coord_collinear()
paired_detail <- magnify_ggchord(
  source + coord_collinear(), regions, detail, side = "bottom"
)
paired_detail
export_ggchord_fusion(paired_detail, width = 12, height = 8)$focus_regions
```

This does not create multiple insets or infer base-level alignments. Each interval must be fully visible in both panels. Instance selectors apply per row.
