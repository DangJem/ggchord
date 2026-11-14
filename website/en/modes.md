## Four modules and their responsibilities {#modules}

| Module | Purpose | Coordinate or composition |
| --- | --- | --- |
| [I · Multi-sequence relationships](multi-sequence.html) | Compare sequence intervals, annotations and links | Chord, Collinear, Radial and Sankey |
| [II · Single-sequence circular maps](circular.html) | Annotate one complete sequence | Circular |
| [III · HSP Dotplot](dotplot.html) | Display query/subject HSP endpoints on two axes | Dotplot |
| [IV · Multi-view composition](fusion.html) | Place independent registered views and connect source intervals | Fusion composition layer |

There are six Coord types. Fusion combines plots; it is not a seventh coordinate system. Dotplot uses a dedicated segment layer and is not a registered Fusion endpoint.

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

## Choose a relationship view {#choose}

### Chord {#chord}
An overview of relationships around curved backbones. [Geometry and curvature](chord-geometry.html).

### Collinear {#collinear}
Aligned rows for sequence order, local windows and anchor comparisons. [Alignment and arrangement](collinear-layout.html).

### Radial {#radial}
Radial spokes retain sequence coordinates while showing syntenic links around the centre. [Radial settings](radial-synteny.html).

### Sankey {#sankey}
Two columns repeat each source as a query and subject display instance. The bands represent observed intervals, not conserved flow. [Sankey settings](sankey.html).

### Circular {#circular}
One sequence ring with prepared annotations, primers and restriction sites. [Circular maps](circular.html).

### Dotplot and Fusion {#composition}
Dotplot places HSP endpoints on two sequence axes. Fusion lays out already-valid independent plots and connects explicit source intervals. Their data and geometry contracts differ.

## Source identity and units {#terms}
Biological strand, display orientation and alignment direction are distinct. bp and aa are separate coordinate domains; cross-domain projection needs an explicit CDS map. Source lengths are never inferred from the last visible feature.

{{terms}}


<span id="check-result"></span><span id="coordinate-controls"></span>

## Get started {#getting-started}

[Install and draw a first plot](installation.html#first-plot), then choose inputs and layers from the four modules above. With existing tables, check [minimal fields](data.html#tables), identifiers, complete lengths and units.

## Shared techniques {#shared-techniques}

[Flexible parameter formats](parameters.html) distinguish shared, named, ordered and strand settings for supported geometry and Position parameters. Use aes/scales for data differences and themes for appearance; see the grammar table above.

[Feature shapes](feature-geometry.html) · [Labels](labels.html) · [Tracks](tracks.html) · [Styling](styling.html) · [Output](output.html) · [Diagnostics](diagnostics.html)
