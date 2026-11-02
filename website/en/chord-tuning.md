These controls have different owners: geoms configure sequence geometry, Positions place local features, coords select the global view, and scales/guides control encodings.

## Inspect directional evidence first {#prepare}
Extract the [input bundle](../assets/tutorial-inputs.zip) first.
{{prepare:sequencepalette,docstyle,showcasedata}}

## Inspect directional evidence {#orientation}
`show_optimization` optimizes display orientation only, preserving sequence order and weighting observed alignment lengths. Inspect seq_orientation and score_before/score_after. A display flip is not a new biological inversion result; the heuristic does not guarantee a global optimum.

## Combine capabilities into a readable view {#parameters}
{{code:showcasechord}}
{{figure:showcasechord}}
Named vectors match sequence IDs. Curvature 0 is straight, signed values bend in opposite directions, and magnitude changes bow depth. Radius changes sequence placement. The result combines automatic orientation, varying radius, zero/positive/negative curvature and local clearance without changing source coordinates.

## Why clearance modes look different {#clearance}
`link_avoid = "uniform"` clears the whole endpoint interval. `"smooth"` retreats locally near feature polygons; `"none"` disables feature avoidance. Automatic clearance requires `link_gap = NULL`; a fixed number overrides it. Clearance considers feature geometry, not text, and does not automatically route around intervening sequence rows.

## Guides and export {#guides}
Identity uses a long left-side colourbar; sequence and strand keys remain below. Explicit strand overrides mirror the arrow glyphs. See [styling](styling.html) for guide responsibilities. Save at 12×10 in, then review [1.7 The complete showcase](chord.html).
