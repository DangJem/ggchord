Adjacent and cross-row links impose different reading costs. Separate their data and controls while preserving observed HSP endpoints. Control points change display paths only.

## Prepare this page {#prepare}
Extract the [input bundle](../assets/tutorial-inputs.zip) first.
{{prepare:sequencepalette,docstyle,showcasedata}}

## Specify controls for long links {#routes}
{{code:showcasecollinear}}
{{figure:showcasecollinear}}
The linear showcase omits the crowded tail assembly protein and gp4 protein candidates, retaining them in the source and chord view; the remaining 16 features are individually labelled. Keep the selected adjacent HSPs and one longest HSP per non-neighbour pair. `show_cross_controls` matches `show_cross` row order; rebuild both after sorting/filtering. Controls use exported layout x/y units, not bp. Each link accepts two or four numbers, shaping its two ribbon edges rather than invoking a global route planner.

Shorter sequences are horizontally centred. A backbone-only layout provides actual row positions, lengths and offsets; controls are placed just beyond the intervening right ends rather than far from the data. Row spacing uses seq_gap and the longest backbone width, while controls do not accept bp coordinates. Changing spacing, order or orientation requires inspecting the paths again; arbitrary controls are not guaranteed free of self-intersection.

## Local clearance and dashed outlines {#clearance}
`position_strand()` separates strands. `link_avoid = "smooth"` with the default NULL gap retreats locally near arrows. Non-adjacent outlines are dashed by default; a larger linewidth makes them visible without overriding linetype. Adjacency follows display order, not biological alignment direction.

## Why the ruler was offset {#ruler}
`geom_seq_scale()` defaults to a 2 mm horizontal and 6 mm vertical offset. The earlier inset was configuration, not incorrect bp conversion. `offset = unit(c(0, 10), "mm")` aligns the start with the left sequence envelope and reserves vertical label space. With centred rows, source coordinate0 is no longer aligned across sequences; the ruler uses the left envelope, and reversed rows have their source origin at the other end. A scale bar measures distance rather than a particular source interval.

## Continue in depth {#next}
Use [focus/segments](views.html) for local questions and [CDS neighbourhoods](homology.html) for homology evidence; review the [2.6 Showcase](collinear.html). GC remains available in [tracks](tracks.html), while the flagship reserves space for labels and links.
