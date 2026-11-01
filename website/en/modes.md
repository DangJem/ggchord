<div class="task-intro"><p>Choose a layout for your question. The examples below compare three views using source-backed bundled records.</p></div>

## Which view should you choose? {#choose}
| Task | Preferred view | Avoid this misreading |
| --- | --- | --- |
| Overall relationships among several sequences | Chord | Crossings alone are not evidence of inversion |
| Annotations on one plasmid or circular sequence | Circular | Multiple sequences do not fit a single-sequence coordinate |
| Positions, gene order and local neighbourhoods | Collinear | Display alignment or flipping is not a biological conclusion |

All three consume existing data; none computes a new BLAST alignment. Try the minimal combinations below, then follow the [raw-record workflow](real-data.html), [complete plasmid map](plasmids.html) or [local homolog comparison](homology.html).

## Prepare this page {#prepare}
In a fresh R session, run the code below, then follow this page in order. Colours and layout are optional display choices; inputs select only the fields needed.

{{prepare:setup,sequencepalette,featurepalette,docstyle,plasmiddata}}

Run this page’s preparation code below. once. Each object below is a normal
composable ggplot; print it or use `ggsave()`.

## Chord: compare multiple sequences {#chord}
{{code:chord}}
{{figure:chord}}

Links are coloured by query sequence ID using the backbone category palette; this does not encode identity or homolog groups. See the [raw-record workflow](real-data.html#chord) for colour mapped to reported identity.
`coord_chord()` controls global rotation and fitting. Sequence-local radius,
curvature, gap and orientation belong to `geom_seq()`. A NULL link control point
uses the historical origin-based chord route; it does not mean a universal
straight link. Ribbon avoidance defaults to `uniform`, line avoidance to `none`.

## Circular: annotate one sequence {#circular}
{{code:circular}}
{{figure:circular}}
This example selects CDS and the replication origin for display, retaining the full sequence bp axis.
`coord_circular()` requires a single sequence. General features can carry unknown
strand. `geom_feature_plasmid()` is a convenience preset with its own defaults;
choose it explicitly when those defaults fit the task. Primer and restriction
layers are described in [recognition](recognition.html).

## Collinear: align local comparisons {#collinear}
{{code:collinear}}
{{figure:collinear}}
Each sequence is anchored to its own longest known gene for display; this does not assert homology.
`anchor` aligns specified intervals; `orientation = "auto"` chooses supported
orientations using available evidence. Inspect orientation diagnostics rather
than assuming biological inversion. Ribbons connect adjacent displayed rows by
default; use `adjacent_only = FALSE` or explicit `pairs` to retain others.
Non-adjacent outlines default to dashed; explicit styles take priority. These
rules differ from chord mode. [Focus and segments](views.html) preserve source
coordinates for cropped or locally reversed views.

## Check your result {#check-result}
First count sequences: use circular for one annotated sequence, chord for the overall link pattern, and collinear for positions, distances or neighbourhoods.

<div class="next-step">Continue: <a href="installation.html">Get started: your first plot →</a></div>
