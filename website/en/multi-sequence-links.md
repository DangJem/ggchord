{{code:courseblockag}}

{{prepare:frontbackbone,frontlinks}}

## Change where a ribbon passes {#control}

{{code:frontbackbone,frontlinks,frontribbon}}

{{figure:frontribbon}}

The default control is at the layout origin. Set `link_ctrl_point` to `c(0, 0.8)` to move the control upward in layout coordinates. Compare the path through the middle. Both ribbon ends still mark the same intervals.

{{code:frontroute}}

{{figure:frontroute}}

Control points use unrotated layout coordinates. With the default rotation of 0°, (0, 0.8) pulls the ribbon upwards. Rotating the whole layout also rotates the control point. A Bézier control pulls the curve without generally lying on it.

## Compare clearance around a real CDS {#clearance}

Take the CDS at 48,053–49,111 bp on `MT108731.1`. It falls within the alignment's query interval. Move the plus-strand arrow inward by `0.05` layout units so it sits on the ribbon side of the backbone.

{{code:frontobstacle,frontnoavoid}}

{{figure:frontnoavoid}}

`link_avoid = "none"` uses fixed endpoint clearance. Look near 49 kb on the `MT108731.1` backbone, where the ribbon overlaps part of the arrow.

{{code:frontsmooth}}

{{figure:frontsmooth}}

`link_avoid = "smooth"` increases clearance locally near the CDS. At the same spot on the `MT108731.1` backbone, the ribbon edge recedes to leave space for the arrow. The default `"uniform"` applies one clearance across each endpoint interval. A numeric `link_gap` sets the spacing directly.

## Link two positions and highlight an interval {#point-markers}

{{code:frontmarkers}}

{{figure:frontmarkers}}

The outer band marks the full interval on the query sequence. `region_side = "outside"` places it outside the backbone. The straight line joins the two interval midpoints, and the `arrow` points to the subject. Use a ribbon to show the full span on both sides.

Next, [place sequence and gene labels](multi-sequence-labels.html#sequence-labels).

## Branches, filtering and fixed-style boundaries {#branch-contract}

| Parameter | Current default | When to change and constraints |
| --- | --- | --- |
| `link_branch` | `"none"` | query/subject share trunks for equal endpoints; unequal endpoints do not establish homology |
| `link_branch_fraction` | `0.15` | Trunk length relative to the shortest member, in (0, 0.5]; only meaningful with branching |
| `pairs` | `NULL` | Explicit sequence-pair selection; does not create alignments |
| `link_gap` | `NULL` | Uses avoidance policy; an explicit number overrides automatic clearance |
| `link_avoid` | `"uniform"` | none uses fixed clearance, smooth local avoidance; text and leaders are ignored |
| `link_ctrl_point` | `NULL` | Origin in chord mode; endpoint-derived in collinear mode; not bp coordinates |
| `nonadjacent_linewidth` | `0.4` | Only nonadjacent collinear borders; 0 reuses ordinary linewidth |

Set fixed `fill/alpha/colour/linewidth` outside `aes()`. Map data differences to dedicated ribbon roles with matching scales. Explicit fixed values override matching mappings; omit fixed alpha when mapping opacity. Colour differences do not establish biological relationships. See [`geom_link_ribbon()`](reference/geom_link_ribbon.html) for the full signature.

Continue with [shared branches and point links](multi-sequence-links.html#input) to compare trunks, trunk fractions and curved/straight paths.

## Cross-row links {#cross-row}

{{prepare:frontlinearbase,frontlinks,frontlinearlink}}

## Allow a cross-row link {#routes}

{{code:frontlinearcross}}

{{figure:frontlinearcross}}

Only Collinear uses thin solid outlines for adjacent rows and heavier dashed outlines for cross-row ribbons. Chord, Radial and Sankey use ordinary solid outlines by default. `linewidth` controls the ordinary width; `nonadjacent_linewidth` separately controls cross-row width and defaults to 0.4.

The second ribbon layer receives one alignment between the first and third rows. Set `adjacent_only = FALSE` to show it. The default dashed outline marks the skipped row. The original ribbon between rows one and two stays in its own layer.

## Set the cross-row outline width {#width}

{{code:frontcrosswidth}}
{{figure:frontcrosswidth}}

`nonadjacent_linewidth = 0.8` thickens the second layer's cross-row ribbon. The first layer's adjacent ribbon and both line types are retained.

{{code:frontcrosssame}}
{{figure:frontcrosssame}}

Zero uses the layer's ordinary `linewidth`, here 0.2. The outline remains dashed. When using `aes(linewidth = ...)`, the mapping and scale control the widths.

## Shared branches and point links {#branches-and-points}

## Prepare an illustrative shared endpoint {#input}

{{code:comparebranchdata}}

Three illustrative sequences each span 1000 bp. Two links share query interval 200–600 bp and connect to different target intervals. Coordinates are one-based inclusive. Every figure uses 9 × 7 inches; endpoints stay fixed across branch settings.

## Independent ribbons and a shared query trunk {#branches}

{{code:comparebranchnone}}
{{figure:comparebranchnone}}

With default `link_branch = "none"`, the ribbons travel independently and may overlap near the query end.

{{code:comparebranchquery}}
{{figure:comparebranchquery}}

Only change to `"query"`: equal query endpoints can share a trunk. Overlapping intervals are not necessarily equal endpoints. Real data without suitable groups may show no change when branching is enabled. `"subject"` groups equal subject endpoints instead.

## Change the trunk fraction {#fraction}

{{code:comparebranchfraction}}
{{figure:comparebranchfraction}}

Relative to the query-branch figure, only `link_branch_fraction` changes from default 0.15 to 0.35. The shared trunk becomes longer relative to the shortest member. It is not identity, a probability or a newly measured alignment length.

| Parameter | Default | Effect and constraints |
| --- | --- | --- |
| `link_branch` | `"none"` | none/query/subject; selects the shared endpoint |
| `link_branch_fraction` | `0.15` | In (0, 0.5]; none does not create a shared trunk |
| `link_ctrl_point` | `NULL` | Shapes far ends when branching; the trunk determines its junction |
| `pairs` | `NULL` | qaccver/saccver table selecting existing pairs, not generating alignments |

## Replace intervals with points {#point-links}

Point links require explicit `qpos/spos`; an interval table is not a point table. Here choose interval midpoints. This loses interval length encoded by ribbon width, so explain the selection rule first.

{{code:comparepointdata}}
{{code:comparelinkcurve}}
{{figure:comparelinkcurve}}

Default `link_type = "curve"` connects the selected positions with curves.

{{code:comparelinkstraight}}
{{figure:comparelinkstraight}}

Only change to `"straight"`, keeping input, anchors and dimensions fixed. Add `arrow = grid::arrow()` separately for direction. It represents the query-to-subject drawing convention, not independent evidence of gene strand or transcription.

## Real inputs and coordinate choices {#adapt}

Preserve real query/subject endpoints. Do not edit unequal endpoints to manufacture shared trunks. See [link preparation](ribbons.html) for filtering and bundling; geometric branching does not perform those statistics.

Chord and collinear modes accept these link layers; single circular mode rejects ribbon/link. Linear adjacency filtering still applies. Check `adjacent_only` or explicit `pairs` before showing cross-row ribbons; see [cross-row links](multi-sequence-links.html#routes).

Example: change only the trunk fraction to 0.25. Return to the default fraction and compare curved versus straight point links. With your own data, first decide whether interval lengths or individual positions matter, then consult [`geom_link_ribbon()`](reference/geom_link_ribbon.html) and [`geom_link_line()`](reference/geom_link_line.html).

{{prepare:ribboninputs}}

## Choose interval ribbons or point links {#interval-points}

{{code:pointlinks}}

{{figure:pointlinks}}

These arrows join the midpoints of the first three alignments and point toward the subject. Use ribbons to show complete intervals.

