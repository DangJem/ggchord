<a id="modules"></a>

## Four modules and their responsibilities {#modules}

| Module | Purpose | Coordinate or composition |
| --- | --- | --- |
| [I · Multi-sequence relationships](multi-sequence.md) | Compare sequence intervals, annotations and links | Chord, Collinear, Radial and Sankey |
| [II · Single-sequence circular maps](circular.md) | Annotate one complete sequence | Circular |
| [III · Sequence alignment Dotplot](dotplot.md) | Display query/subject HSP endpoints on two axes | Dotplot |
| [IV · Multi-view composition](fusion.md) | Place independent registered views and connect source intervals | Fusion composition layer |

There are six Coord types. Fusion combines plots; it is not a seventh coordinate system. Dotplot uses a dedicated segment layer and is not a registered Fusion endpoint.

<a id="prepare"></a>

## Shared plotting grammar {#prepare}

| Component | Responsibility | Entry |
| --- | --- | --- |
| Data | Source identifiers, complete lengths, intervals and units | `ggchord()` |
| Geometry | Backbones, links, features, tracks and text | `geom_*()` |
| Mapping | Associate data with an aesthetic | `aes()` |
| Scale | Define colours, ranges and guides | `scale_*()` |
| Position | Place features and tracks locally | `position_*()` |
| Coordinate | Arrange the complete view | `coord_*()` |
| Theme | Non-data appearance | `theme_*()` |

<a id="choose"></a>

## Choose a relationship view {#choose}

<a id="chord"></a>

### Chord {#chord}
An overview of relationships around curved backbones. [Geometry and curvature](chord-geometry.md).

<a id="collinear"></a>

### Collinear {#collinear}
Aligned rows for sequence order, local windows and anchor comparisons. [Alignment and arrangement](collinear-layout.md).

<a id="radial"></a>

### Radial {#radial}
Radial spokes retain sequence coordinates while showing syntenic links around the centre. [Radial settings](radial-synteny.md).

<a id="sankey"></a>

### Sankey {#sankey}
Two columns repeat each source as a query and subject display instance. The bands represent observed intervals, not conserved flow. [Sankey settings](sankey.md).

<a id="circular"></a>

### Circular {#circular}
One sequence ring with prepared annotations, primers and restriction sites. [Circular maps](circular.md).

<a id="composition"></a>

### Dotplot and Fusion {#composition}
Dotplot places HSP endpoints on two sequence axes. Fusion lays out already-valid independent plots and connects explicit source intervals. Their data and geometry contracts differ.

<a id="terms"></a>

## Source identity and units {#terms}
Biological strand, display orientation and alignment direction are distinct. bp and aa are separate coordinate domains; cross-domain projection needs an explicit CDS map. Source lengths are never inferred from the last visible feature.

{{terms}}


<span id="check-result"></span><span id="coordinate-controls"></span>

<a id="getting-started"></a>

## Get started {#getting-started}

[Install and draw a first plot](installation.md#first-plot), then choose inputs and layers from the four modules above. With existing tables, check [minimal fields](data.md#tables), identifiers, complete lengths and units.

<a id="shared-techniques"></a>

## Shared techniques {#shared-techniques}

[Flexible parameter formats](parameters.md) distinguish shared, named, ordered and strand settings for supported geometry and Position parameters. Use aes/scales for data differences and themes for appearance; see the grammar table above.

[Feature shapes](feature-geometry.md) · [Labels](labels.md) · [Tracks](tracks.md) · [Styling](styling.md) · [Output](output.md) · [Diagnostics](diagnostics.md)
