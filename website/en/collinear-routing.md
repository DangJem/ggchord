Collinear plots show links between neighbouring rows by default. Use one observed alignment from the first to the third row to see how to display a cross-row ribbon.

## Start with an adjacent link {#prepare}

{{code:frontlinearbase,frontlinks,frontlinearlink}}

{{figure:frontlinearlink}}

## Allow a cross-row link {#routes}

{{code:frontlinearcross}}

{{figure:frontlinearcross}}

Adjacent ribbons use thin solid outlines and non-adjacent rows use heavier dashed outlines. `linewidth` controls the ordinary width; `nonadjacent_linewidth` separately controls cross-row width and defaults to 0.4.

The second ribbon layer receives one alignment between the first and third rows. Set `adjacent_only = FALSE` to show it. The default dashed outline marks the skipped row. The original ribbon between rows one and two stays in its own layer.

## Set the cross-row outline width {#width}

{{code:frontcrosswidth}}
{{figure:frontcrosswidth}}

`nonadjacent_linewidth = 0.8` thickens the second layer's cross-row ribbon. The first layer's adjacent ribbon and both line types are retained.

{{code:frontcrosssame}}
{{figure:frontcrosssame}}

Zero uses the layer's ordinary `linewidth`, here 0.2. The outline remains dashed. When using `aes(linewidth = ...)`, the mapping and scale control the widths.

## Control points and clearance {#clearance}

If a cross-row ribbon covers the middle row, use `link_ctrl_point` to lead it to the side. Controls use layout x/y coordinates and the alignment endpoints still use bp. See the [four-sequence case study](collinear.html#compose) for per-ribbon controls combined with `link_avoid = "smooth"`.

## Scale bars and coordinate zero {#ruler}
The default bar has a 0 mm horizontal inset and 6 mm vertical spacing. Its start follows the left sequence envelope. Forward, left-aligned rows share that zero point. Read individual zero ticks after centring or reversing rows.

## Continue in depth {#next}
Use [focus/segments](views.html) for local questions and [CDS neighbourhoods](homology.html) for homology evidence; review the [Showcase](collinear.html). GC remains available in [tracks](tracks.html), while the flagship reserves space for labels and links.
