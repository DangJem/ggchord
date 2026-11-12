## Connect real source intervals {#intervals}

Fusion combines independent genomic plots and connects their real sequence
intervals. The center and peripherals can use registered genomic Coord views; Sankey and Radial retain their explicit display instances. Nested Fusions and registered character details are supported. Each plot remains usable on its own. The composition
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
Dotplot, facets, external ggplot coordinate guessing and automatic global routing
remain outside this interface. Explicit derived views and placement are described below. Links may cross or pass through panel interiors; foreground
annotations are retained, but arbitrary layouts are not guaranteed collision-free.


## A source interval and its local detail {#magnify}

`magnify_ggchord()` composes a Circular/Chord/Collinear source and an explicit Collinear or registered character detail plot. A source interval is biological data, not a rectangle in image pixels. The same region must be fully visible in both panels. Both rebuild at the final device size; no screenshot is stretched.

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

The detail plot owns its data, scales, layers and statistics. For count/bin/coverage tracks, prepare and state the intended local or global statistical population yourself; the wrapper records `explicit_detail_plot` and never silently chooses a scope. It does not copy arbitrary source layers, infer alignments or accept unregistered Cartesian character panels. Use one ordered region, explicitly select `instance`/`detail_instance` when an accession is repeated, and declare `wrap = TRUE` for origin-crossing selections as described below. Explicit nested/free-position allocations are described below; measured obstacle routing is described below.

## Paired source windows with real HSPs {#paired-detail}

In this example, each row of `region` selects one original source interval; all rows share one explicit detail panel. The example selects the two endpoints of a real HSP from the four-sequence dataset. Gene direction and the local HSP come from the focused input tables. `region_row` links each derived region, front and decorative leader to its input row; arbitrary original metadata stays in `focus_region_data`.

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

## Independent detail panels {#multiple-details}

Split the two endpoints of one real HSP into independent local views. A named
list supplies the detail plots, `detail_panel` routes each source row, and a
named `side` vector assigns distinct sides. Each panel owns its data and scales.

```r
library(ggchord)
library(dplyr)
seqs <- seq_data_example |> select(accver, length)
hsp <- ribbon_data_example |> slice(1)
regions <- bind_rows(
  transmute(hsp, accver = qaccver, start = pmin(qstart, qend),
    end = pmax(qstart, qend), detail_panel = "query"),
  transmute(hsp, accver = saccver, start = pmin(sstart, send),
    end = pmax(sstart, send), detail_panel = "subject")
)
make_detail <- function(region) {
  local <- focus_ggchord_data(seqs, loci = region)
  ggchord(local$seq_data) +
    geom_seq(show.legend = FALSE) +
    coord_collinear(seq_axis = "vertical")
}
source <- ggchord(seqs, hsp) +
  geom_seq(show.legend = FALSE) +
  geom_link_ribbon(show.legend = FALSE)
details <- list(
  query = make_detail(regions[1, ]),
  subject = make_detail(regions[2, ])
)
multiple <- magnify_ggchord(
  source, regions, details,
  side = c(query = "left", subject = "right")
)
export_ggchord_fusion(multiple, width = 12, height = 8)$focus_regions
```

This example uses four named details, one per cardinal side; diagonal sides and explicit boxes are available below. Every panel needs a region; missing
or unknown targets fail. Source and detail instance selectors remain independent.
`detail_size` accepts a named list of width/height fractions; NULL uses Fusion
defaults. Overlaps fail rather than displacing explicit allocations. Standard
saving rebuilds panels and leaders at the requested size. At small sizes, inspect
labels and enlarge the canvas as needed. Leaders may cross the main plot or one
another; global avoidance is separate from the explicit nested/box allocation described below.

## Character details with exact source anchors {#character-detail}

A Cartesian character plot can be an independent detail panel after `register_ggchord_alignment()`. Register the actual displayed character rows: x is the original alignment/coding column, y is the alignment row ID, and label is the actual symbol. The source map remains separate from the ordinary ggplot scales. These real 60-column excerpts retain the upstream gaps and mismatches; their macro lengths refer to the **ungapped aligned fragments**, not complete genomes.

```r
library(dplyr)
cells <- read_aligned_fasta("examples/alignment/opuntia-aligned.fasta", "DNA")
rows <- unique(cells$accver)[c(1, 3)]
window <- cells |>
  filter(accver %in% rows, between(alignment_column, 321, 380))
sequences <- cells |>
  filter(accver %in% rows, !is_gap) |>
  group_by(accver) |>
  summarise(length = max(source_position), .groups = "drop")
regions <- window |>
  filter(!is_gap) |>
  group_by(source_accver) |>
  summarise(start = min(source_position), end = max(source_position), .groups = "drop") |>
  rename(accver = source_accver)
letters <- ggplot2::ggplot(window, aes(alignment_column, accver, label = symbol)) +
  geom_residue(aes(fill = symbol), size = 2.7, overflow = "error") +
  ggplot2::scale_fill_hue(l = 85, c = 40, guide = "none") +
  ggplot2::theme_minimal()
letters <- register_ggchord_alignment(letters, window)
macro <- ggchord(sequences) + geom_seq() + coord_collinear()
character_detail <- magnify_ggchord(
  macro, regions, letters,
  center_size = c(.22, .7), detail_size = c(.5, .8)
)
export_ggchord_fusion(character_detail, width = 16, height = 5)$focus_regions
# ggsave("character-detail.png", character_detail, width = 16, height = 5)
```

The ordinary plot remains editable: local colours, layers and annotations remain independent. Identity/reverse x scales and discrete y rows are supported. Each Fusion build verifies the registered residue positions/labels after edits and remeasures the final panel in mm; half-cell boundaries attach to the exact endpoint cells. Gap cells have no source position and are never interpolated into bases. A source interval must contain every selected base in biological order; split introns into separate source rows; origin joins use the explicit wrapping declaration below. Use `detail_instance` to select an alignment row when several rows share a source. Unregistered ggplots, facets, nonlinear/flipped coordinates and moved registered cells fail.

`translate_ggchord_cds()` nucleotide rows can be registered in the same way. Keep the translation on codon centres as in `examples/m2-detail.R`; use `source_positions` for full codon membership across segments. AA residue rows use their own aa coordinate domain; connecting a genomic bp macro directly to aa positions still requires an explicit CDS projection. Decorative magnify leaders remain separate from biological Fusion links. Character endpoints have no backbone colour for a Fusion gradient: supply explicit cross-panel fill. Narrow panels retain `geom_residue()`'s tile/error policy; enlarge the physical export instead of relying on DPI. Global leader avoidance is separate; explicit nested/free-position allocation is supported below.

## Eight directions and independent physical placement {#advanced-placement}

`fusion_panel()` supports the four cardinal sides plus `topright`, `topleft`, `bottomright` and `bottomleft`. Default corner layouts reserve space for all eight; side names remain unique. For other arrangements use `angle`/`distance` or an explicit `box`. There is no fixed panel-count limit for explicit boxes/angles: every allocation must still fit. Placement changes the container, not the child's Coord or text orientation. Angles describe physical directions on rectangular exports; distance is a fraction of half the shorter allocation dimension, measured from the centre allocation.

```r
data(seq_data_example)
sequences <- seq_data_example |>
  dplyr::select(accver, length)
base_panels <- ggchord(sequences) + geom_seq(show.legend = FALSE)
central <- base_panels + coord_chord()
mirrored <- base_panels + coord_sankey()
free_angle <- ggchord_fusion(central) +
  fusion_panel(mirrored, "mirror", angle = 31, distance = .8, size = c(.16, .22))
free_box <- ggchord_fusion(central, center_box = c(.02, .35, .2, .8)) +
  fusion_panel(mirrored, "mirror", box = c(.45, .98, .1, .9))
physical_box <- ggchord_fusion(central, center_box = unit(c(5, 85, 10, 100), "mm")) +
  fusion_panel(mirrored, "mirror", box = unit(c(100, 210, 10, 100), "mm"))
export_ggchord_fusion(physical_box, width = 12, height = 6)$panels
```

A numeric box uses xmin, xmax, ymin, ymax fractions; an absolute-unit box keeps its physical dimensions. Coordinates refer to the composition area excluding its cross-panel guide. Do not combine box with side/size/angle/distance, or center_box with center_size. Physical allocations are checked at every target size; they are never silently clamped. Explicit `allow_overlap = TRUE` on `ggchord_fusion()` or `magnify_ggchord()` permits deliberate overlays. Later panels cover earlier panels, using vector grid clipping rather than PNG crops; translucent connections are not repainted twice in the same area. The export's `overlaps` table records covering order. This permits intentional masks and does not promise automatic obstacle avoidance.

Sankey source IDs occur in both mirrored columns. Choose an exact `q_instance`/`s_instance` from `get_chord_layout(plot)$display_instances`; omission fails when ambiguous. Radial and Chord use their own registered backbone frames. All q/s intervals, units and source records remain unchanged. If a declared `source_hsp_id` appears repeatedly, every participating Fusion layer must opt into `allow_repeated = TRUE`, including layers inside nested compositions. `source_observations` distinguishes explicit evidence IDs from display counts. This does not infer identity by comparing equal table values. Dotplot has independent axes and is not yet a registered Fusion endpoint.

## Nested source regions retain their own namespaces {#nested-detail}

A detail may itself be a Fusion. This example uses two explicit windows on the unchanged pUC19 source. The inner source is the larger window; its child is the smaller window. The outer selection targets the inner source through `detail/center`, not an arbitrary rectangle in an image.

```r
data(plasmid_example_pUC19)
source_genome <- plasmid_example_pUC19 |>
  dplyr::select(accver, length)
large_region <- data.frame(accver = "pUC19", start = 101, end = 1900)
small_region <- data.frame(accver = "pUC19", start = 301, end = 600)
large_window <- focus_ggchord_data(source_genome, loci = large_region)
small_window <- focus_ggchord_data(source_genome, loci = small_region)
large_plot <- ggchord(large_window$seq_data) +
  geom_seq(show.legend = FALSE) + coord_collinear()
small_plot <- ggchord(small_window$seq_data) +
  geom_seq(show.legend = FALSE) + coord_collinear()
inner <- magnify_ggchord(large_plot, small_region, small_plot)
whole <- ggchord(source_genome) + geom_seq() + coord_circular()
large_region$detail_panel <- "detail/center"
outer <- magnify_ggchord(
  whole, large_region, inner,
  center_box = c(.02, .22, .15, .85), detail_box = c(.27, .99, .1, .9)
)
report <- export_ggchord_fusion(outer, width = 14, height = 7)
report$focus_regions
report$focus_region_data_by_composition
```

`composition_id` plus `region_row` connects derived regions/fronts/leaders to each original selection in `focus_region_data_by_composition`. Top-level `focus_region_data` remains the original outer input. Nested leaf names use paths such as `detail/center` and `detail/detail`; a composition container is never treated as one genomic axis. If an accession occurs in multiple leaves, select the full panel path and exact instance. `source_panel` can select a source leaf under `center`. Layout exports include actual mm boundaries and requested boxes/angles.

Each plot keeps its own layers, scales and statistical population; the compositor does not silently copy layers or recompute statistics. Origin wraps use the explicit source-length declaration below. The complete negative-strand CDS-to-genome-focus-to-real-character example is in `examples/m3-nested.R`. Registered characters retain their physical tile/error policy; a small nested allocation cannot make 60 letters readable merely by increasing DPI. Full multi-device publication acceptance and broader layout optimisation remain separate concerns.

## Explicit origin-crossing selection {#origin-detail}

Circular layout is not evidence of a circular source. Set `wrap = TRUE` on a
selection with `start > end`. Both selected leaf panels must retain the same
full source length. Character adapters additionally require `circular = TRUE`
and `source_data` with a full sequence or an exact `source_length`; without the
full string, characters remain unverified. At most one source tour is allowed,
so inverse source-to-column mapping stays unique. Negative strands wrap in
their own direction and retain gaps. Different focused instances are never
joined implicitly.

Run the complete real pUC19 workflow in `examples/origin-detail.R`. It extracts
82 observed bases (2640–2686, then 1–35) from the unchanged bundled FASTA.
`focus_region_data` retains the original single row. Derived tables identify
two pieces with `region_piece`, sharing one `region_row`; only the two external
boundaries receive leaders. The seam is not another selection or relationship.

```r
# Requires dplyr and the source checkout or tutorial input bundle.
source("examples/origin-detail.R")
origin_detail
```

## Measured routes and display corridors {#routing}

`fusion_ribbon(route = "avoid")` keeps both complete interval fronts and their
physical departure normals. It routes a narrow polygonal body around measured
other leaf containers, native layer text, and final gtable guide/decoration
cells. `route_width` is in mm and does not encode an alignment length or count.
The default `"direct"` preserves the existing Bezier drawing and reports its
collisions; it does not silently filter records.

Use `route = "bundle"` with an explicit character `route_group` column to share
a corridor within one layer and ordered panel pair. Each source row still has
its own front, colour/alpha and polygon. Shared-channel overpainting is not
read support or another observation. Groups are user declarations; no cluster
or relationship is inferred. Run the real three-HSP example:

```r
# Requires dplyr and the source checkout or tutorial input bundle.
source("examples/fusion-routing.R")
routed
```

`magnify_ggchord(route = "avoid")` applies the same measured-space contract to
its decorative leaders. Exported `leader_paths` contains the final polyline
vertices, separate from biological `paths`. `routing_corridors` identifies
shared channels and member source rows; `routing_obstacles` records the
physical envelopes. Connections/leaders report route status and collisions.
Each target size is rebuilt, including ancestor checks after nested allocation.

No clear route produces an explicit error: move/enlarge the panels, narrow a
selection, revise a declared group, or inspect `route = "direct"`. Ancestor
overlays that obstruct an already measured child route also fail explicitly.
Routing uses deterministic visibility paths, not a global layout optimum.
Rotated text uses conservative bounding rectangles, curved strings larger
envelopes, and native decoration cells their final gtable bounds. Arbitrary
custom grobs and every endpoint data shape are not claimed to be glyph-level
obstacles. Theme changes affect metrics, never source coordinates or evidence.

## Outer Guides and the panel budget {#outer-guides}

Keep each panel's local scales and Guides independent. Composition scales and
`guides()` describe the cross-panel ribbons. Native Guides may occupy different
outside positions, rather than all requiring one shared side:

```r
library(ggchord)
guide_source <- ggchord(seq_data_example) + geom_seq(show.legend = FALSE)
guide_fusion <- ggchord_fusion(guide_source, layout = "adaptive") +
  fusion_panel(guide_source + coord_collinear(), "detail") +
  fusion_ribbon(ribbon_data_example, aes(fill = qaccver, alpha = pident),
    qpanel = "center", spanel = "detail") +
  ggplot2::guides(
    fill = ggplot2::guide_legend(position = "bottom", nrow = 1),
    alpha = ggplot2::guide_legend(position = "right")
  )
guide_fusion
```

`layout = "fixed"` keeps existing placements. `"adaptive"` reuses the native
Guide candidate/measurement engine, jointly considers side budgets and panel
allocations, and then rebuilds every child in the remaining physical space.
It may move/reflow unlocked Guides; explicit position, direction, rows/columns,
fonts, key dimensions and absolute panel boxes remain constraints. For example,
a default Identity ruler can become horizontal on a wide, short canvas; its
length is remeasured for all six original labels rather than dropping 95 or 100.

`metadata$guide_layout` reports choices and locks. `guide_bounds_mm` and
`guide_sizes_mm` describe the final outer boxes, after packaging and before
panel measurement. Nested composition Guide reports appear in `diagnostics`.
Impossible constraints fail explicitly. `theme(legend.position = "none")`
hides composition Guides using native ggplot2 semantics. Inside composition
Guides remain unsupported; each independent panel can still own local Guides.
This is a bounded physical-space search, not a claim of global optimality.
