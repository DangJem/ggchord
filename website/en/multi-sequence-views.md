<a id="comparison"></a>

## Follow one alignment through four coordinate systems {#comparison}

All four plots use the full lengths of MT108731.1 (64,323 bp) and MT118296.1 (32,090 bp) and the first local alignment between them in that direction. Each example runs on its own. The same layers supply references and alignment intervals. Each view chooses its own spacing and coordinate system.

<a id="chord"></a>

### Find the two endpoint intervals on arcs {#chord}

{{example:views-chord}}

The longer blue arc belongs to MT108731.1. The band joins 26,298–51,139 bp on that arc to 7,121–31,959 bp on the shorter orange arc. The arrow on each backbone shows its display direction. [Build this diagram a layer at a time](chord.md).

<a id="collinear"></a>

### Read both intervals against straight axes {#collinear}

{{example:views-collinear}}

The shorter sequence now lies below the longer one. Both rows begin at base 1 on the left. The lower end of the band reaches nearly to the right end of its backbone. [Add the other two alignments](collinear.md#adapt).

<a id="radial"></a>

### Read from the centre out along each spoke {#radial}

{{example:views-radial}}

Here base 1 on each spoke faces the central gap. Follow the band between the upper and lower spokes to find the same two intervals. [Add a third spoke](radial-synteny.md#adapt) to see how two alignments can share a query.

<a id="sankey"></a>

### Follow query on the left to subject on the right {#sankey}

{{example:views-sankey}}

Each accession appears in both columns. The single input row links the upper left copy to the lower right copy. Next, [compare the three bands for this pair](sankey.md#adapt), or [use their coordinates on two dotplot axes](dotplot.md).
