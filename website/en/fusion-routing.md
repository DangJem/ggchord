{{code:courseblockbk}}

{{prepare:fusiondata,fusionpanels,fusionlinear,fusionview}}

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

The first example uses one panel per named side: `right`, `top`, `left`, `bottom`. Diagonal, free and automatic placement are described below. IDs must be unique;
`center` is reserved. `size=c(width,height)` and `center_size` are fractions of
the content area after native outer decoration and cross-panel Guides have been allocated. Explicit
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
the source workflow explicitly. Ordinary `fusion_ribbon()` inputs require prepared continuous intervals.
For cross-origin selections, use the explicit `magnify_ggchord()` with region `wrap = TRUE`
region workflow below, with matching complete source lengths. Direct bp↔aa
links, Dotplot, facets and unregistered external plots remain unsupported.
`direct`, `avoid` and `bundle` route against their documented measured obstacles;
they do not guarantee a global optimum. Explicit derived views and placement are described below. Links may cross or pass through panel interiors; foreground
annotations are retained, but arbitrary layouts are not guaranteed collision-free.
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

{{code:courseblockbl}}

![Explicit display corridors and obstacles](../assets/courserouted.png)

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
