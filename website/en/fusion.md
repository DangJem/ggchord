## Connect real source intervals {#intervals}

Fusion combines independent genomic plots and connects their real sequence
intervals. The center can be Circular or Chord; peripherals are horizontal or
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
