First change the path with a control point. Then compare endpoint clearance near a CDS using the same ribbon. These settings act on the route through the middle and the space near the endpoints, respectively.

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

The outer band marks the full interval on the query sequence. `region_side = "outside"` places it outside the backbone. The straight line joins the two interval midpoints, and the arrow points to the subject. Use a ribbon to show the full span on both sides.

Next, [place sequence and gene labels](chord-labels.html).

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

Continue with [shared branches and point links](chord-branches.html) to compare trunks, trunk fractions and curved/straight paths.

## Region side, width and offset {#region-parameters}

Keep the original query ID and interval of the first observed HSP. Defaults are width 0.08, offset 0 and inside. Each subsequent figure changes only outside placement, width or offset. These are local geometry quantities, not bp. Set fixed fill in geom parameters or map categories to region_fill with its dedicated scale. Highlighting does not establish a new functional annotation.

{{code:comparefinishchord}}
{{code:compareregioninside}}
{{figure:compareregioninside}}
{{code:compareregionoutside}}
{{figure:compareregionoutside}}
{{code:compareregionwidth}}
{{figure:compareregionwidth}}
{{code:compareregionoffset}}
{{figure:compareregionoffset}}
