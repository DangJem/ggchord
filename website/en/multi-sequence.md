<a id="task"></a>

## Shared data, four views {#task}

Complete sequence lengths, observed HSPs and real CDS can share four relationship coordinates. The [strict comparison](multi-sequence-views.md#input) keeps inputs, layers and mappings fixed with explicit all-pair selection. The Chord flagship below has its own display choices; see [flagship notes](multi-sequence-views.md#flagship).

{{manualfigure:flagchord}}

<a id="path"></a>

## What do you want to draw or adjust? {#path}

- [Draw backbones, links and annotations](multi-sequence-start.md#draw)
- [Compare four layouts with the same inputs](multi-sequence-views.md#input)
- [Change sequence order and display orientation](multi-sequence-layout.md)
- [Adjust link endpoints and routes](multi-sequence-links.md)
- [Show sequence names and feature labels](multi-sequence-labels.md)

<a id="chord"></a>

## Chord {#chord}
Curved backbones, radius, curvature and routed interval links. [Settings](chord-geometry.md) · [Example](chord.md).

<a id="collinear"></a>

## Collinear {#collinear}
Row or column arrangement, anchors, orientation and local windows. [Settings](collinear-layout.md) · [Example](collinear.md).

<a id="radial"></a>

## Radial {#radial}
Sequence spokes, start angle and angular arrangement. [Settings and source instances](radial-synteny.md).

<a id="sankey"></a>

## Sankey {#sankey}
Mirrored query/subject columns, column gap and repeated display instances. [Settings](sankey.md).

<a id="inputs"></a>

## Input contract {#inputs}
Complete source lengths, identifiers and units remain independent of display geometry. Biological strand comes from annotation; changing the coordinate does not rewrite it. [Data contract](data.md) · [Compatibility](compatibility.md).
