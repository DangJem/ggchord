Use the same four backbones and one alignment to compare row spacing and add a scale bar.

## Default layout {#baseline}

{{code:frontlinearbase,frontlinks,frontlinearlink}}

{{figure:frontlinearlink}}

## Increase row spacing {#gap}

{{code:frontlineargap}}

{{figure:frontlineargap}}

In a collinear plot, `seq_gap` is relative to the longest backbone width and defaults to `0.30`. Set the gap below MT108731.1 to `0.45`. The other gaps retain `0.30`; the final row has no following gap. The ribbon still joins the same bp intervals.

## Add a fixed-length scale bar {#scale}

{{code:frontlinearscale}}

{{figure:frontlinearscale}}

Return to the default spacing and add a 10 kb scale bar. It spans about a third of the 32 kb backbone. Next, see [inputs and display orientation](collinear-data.html).
