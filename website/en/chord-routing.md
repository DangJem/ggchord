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
