## Plot original HSP endpoints {#dotplot}

The source-backed EMALE PAF pair has 11 observed HSPs. Its declared sequence
lengths supply both axes. `read_paf()` converts zero-based half-open intervals
to the internal one-based convention and reverses target endpoints for minus
alignments. The files and their provenance are in the [input download](formats.html).
The following code needs the optional `dplyr` package for data preparation.

{{code:dotplotdata}}

Each line runs from `(qstart, sstart)` to `(qend, send)`. `geom_seq()` is not
required. Use `link_*` mappings and scales for colour, alpha, width and type.
The default Identity domain is 0–100 with ticks 0/50/80/90/95/100; an explicit
user scale overrides it.

{{code:dotplotview}}

{{figure:viewdotplot}}

The legacy constructor may warn about descending endpoints in these valid reverse HSPs. Keep the normalized ordering; sorting each interval would erase the displayed inversion.

`x_seqs` and `y_seqs` independently select/order accessions. `gap` is display
spacing in native bp/aa, not omitted source sequence. Axis directions accept
1/-1 or complete named vectors. Display reversal changes the slope without
changing the biological source table. Each axis starts at source position 1;
concatenated segments retain their own accession and local tick labels.

{{code:dotplotreverse}}

Use `map_ggchord_position(..., axis = "x"/"y")` for native/global conversion,
with `inverse = TRUE` to go back. Display gaps return NA. Focused source
origins use `source = TRUE`; default axis labels stay window-relative.
Unknown directions are warned, omitted and reported; contradictory direction
metadata is rejected. Repeated declared HSP IDs require `allow_repeated=TRUE`
in all participating layers. Unidentified duplicate/self rows remain distinct.
Without cg/cs/CIGAR or another alignment trace, these are HSP segments, not
base-level mismatches or gap paths.

## Proportional radial spokes {#radial}

The bundled four-sequence data retain their original lengths and CDS strands.
This explicit display filter selects observed HSPs of at least 900 bp and
75 percent identity, taking the two longest per directed pair. The count track
counts those selected query-start observations, not inferred depth or coverage.

{{code:radialdata}}

{{code:radialview}}

{{figure:viewradial}}

Spoke order comes from `geom_seq(seq_order=...)`. Default lengths are
proportional; `length_mode="equal"` explicitly normalizes display lengths
while preserving source coordinates. `start_angle` is measured counterclockwise
from positive x; `clockwise` controls spoke order, not source strand. Positive
feature/track offsets point counterclockwise from an outward spoke. Physical
Track widths and source statistics reuse the shared band and Stat system.

Cyclic neighbors include the last/first pair **only when an HSP exists**.
Without that relationship its channel stays empty. N=2 draws each source row
once; N=1 is rejected. Reverse HSPs retain full ordered interval fronts and
twists. New Radial fronts use inclusive interval boundaries `[min-1,max]`, so
a one-position HSP has nonzero width; axis positions stay one-based. Existing
views retain their own documented interval conventions. The hole is display
spacing, not deleted sequence. Sequence names sit beyond the outer tips.

{{code:radialrotate}}

Ribbons can select nonneighbors with explicit `pairs` or
`adjacent_only=FALSE`. Nonneighbor point links need `between`; named A–B–A
instances use the same unambiguous [routing contract](foundation.html#instances).
User control points use final x/y coordinates. Endpoint gradients and subject
tips remain explicit options; default ribbons retain Identity fills.

## Export and reverse-map {#export}

{{code:radialmap}}

Exports retain `radial_frames`, original source intervals and projected complete
fronts in `radial_endpoints`, plus connection selection, instance routing,
track statistics and label visibility. Radial inverse mapping requires both
x and y on the selected backbone. This avoids inventing a position from a
single x value shared by multiple spokes. Focus/segment inputs retain source
origins; radial display does not join multiple windows onto one spoke.

## Supported combinations {#support}

| Layer or operation | Dotplot | Radial |
| --- | --- | --- |
| `geom_alignment_segment()` | Identity HSP segments | Rejected |
| Sequence, gene/feature, region and break | Rejected | Source-local projection |
| Manual/automatic gene/feature labels | Rejected | Shared measured queues; inspect hiding at target size |
| `geom_track()` and existing Track shortcuts | Rejected | Shared source Stat and band projection |
| Ribbon/point links | Rejected | Cyclic/default or explicitly routed relationships |
| Native ggplot2 annotations | Global concatenated units | Final layout x/y units |
| Shared branches, facets | Rejected | Rejected |
| Circular primers/restriction sites/centre labels | Rejected | Rejected |
| Shared straight genomic scale bar | Rejected | Rejected; use each spoke's native axis |

These are supported contracts and bounded examples, not every possible
parameter combination. Central clearance conflicts produce an actionable
error; explicit control points may cross the hole or other spokes. Dense
labels can be hidden by layer policy, with source identities retained in
exports. Complex Fusion and fine alignment traces belong to later stages.
