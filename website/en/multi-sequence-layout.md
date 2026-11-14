{{code:courseblockaf}}

{{prepare:frontbackbone}}

<a id="direction"></a>

## Orientation: start the coordinates at the other end {#direction}

{{code:frontorientation}}

{{figure:frontorientation}}

`seq_orientation = -1` reverses the sequence's display direction. Compared with the default, this backbone’s zero tick and arrow swap ends. Input coordinates and biological strands still come from the tables. The three vectors are named by sequence ID, so they still match if you reorder the sequences with `seq_order`.

<a id="order"></a>

## Order: arrange the backbones {#order}

{{code:frontorder}}

{{figure:frontorder}}

`seq_order` takes the sequence IDs in display order. Here `OQ646790.1` moves ahead of `MT118296.1`. Link rows still match the same IDs after reordering.

Next, [change the path and endpoint clearance of a ribbon](multi-sequence-links.md#control).

<a id="anchors"></a>

## Linear anchors {#anchors}

[Collinear](collinear-layout.md#anchor) · [Radial](radial-synteny.md) · [Sankey](sankey.md)
