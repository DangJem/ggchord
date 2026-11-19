<a id="comparison"></a>

## Keep the observations fixed {#comparison}

All four examples use the same two complete references and the same three observed HSPs from MT108731.1 to MT118296.1. Each code block runs independently. Compare where the ends of each ribbon meet a backbone. The order of HSP endpoints and the lengths of both reference sequences are preserved across these views.

<a id="chord"></a>

### Chord {#chord}

Read ribbon endpoints around the circumference. This layout gives many-to-many relations a compact overview.

{{example:views-chord}}

<a id="collinear"></a>

### Collinear {#collinear}

Read along each row to compare the positions and spans of the three alignments. To change the displayed direction of one reference, see [order and orientation](multi-sequence-layout.md).

{{example:views-collinear}}

<a id="radial"></a>

### Radial synteny {#radial}

Read outward along each spoke. With additional references, the spoke order determines which neighbours sit next to each other.

{{example:views-radial}}

<a id="sankey"></a>

### Sankey {#sankey}

Read the query copy on the left and the subject copy on the right. Ribbon width at each end reflects the observed span on that sequence. The two-column arrangement is useful when query and subject roles matter.

{{example:views-sankey}}

The [gallery](gallery.md) contains full multi-sequence examples with their data and rendering scripts.
