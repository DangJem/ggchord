In radial synteny, each reference is a straight spoke with its own position scale. Follow a ribbon to the two spokes to find the intervals it connects.

<a id="draw"></a>

## Read one HSP on two spokes {#draw}

Use the first alignment from MT108731.1 to MT118296.1. The data and layers match the [two-row comparison](collinear.md). Use the radial projection and its spoke layout for this view.

{{example:views-radial}}

With two references the spokes point in opposite directions. Read outward from the central gap on each one. The subject end of the ribbon reaches almost to the outer end of MT118296.1 at 32,090 bp.

<a id="adapt"></a>

## Add a second subject spoke {#adapt}

Add OQ646790.1 to the references. Group hits by subject before selecting the first row so each subject contributes one observed interval pair. MT108731.1 now connects to two separate spokes, with different widths at the query end.

{{example:first-radial}}

The first spoke starts at noon. In larger maps, set `start_angle` in `coord_radial_synteny()` to rotate the spokes, and `inner_radius` in `layout_sequences()` to change the central opening. Next, [separate query and subject into columns](sankey.md), or [add genes to a four-spoke map](case-radial.md).
