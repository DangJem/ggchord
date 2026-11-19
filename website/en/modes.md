<span id="input"></span><span id="read"></span><span id="save"></span><span id="source"></span><span id="stable"></span><span id="verified-versions"></span>
<a id="modules"></a>

## What do you want to plot? {#modules}

| I want to… | Read the guide | Try a plot |
| --- | --- | --- |
| Compare sequence order, annotations and interval links | [Multi-sequence relationships](multi-sequence.md) | [Start with sequence backbones](multi-sequence-start.md#draw) |
| Show features and tracks along one circular sequence | [Single-sequence circular maps](circular.md) | [Start with lengths and known annotations](circular-start.md#input) |
| Inspect local alignments on two axes | [Alignment dot plots](dotplot.md) | [Draw HSP endpoints and directions](dotplot-start.md#dotplot) |
| Combine an overview, a magnified region or different views | [Multi-view composition](fusion.md) | [Combine and connect two plots](fusion-start.md#panels) |

Browse the [example gallery](gallery.md) by plot type or task for reproducible examples.

<a id="choose"></a>

## Which layout should I use to compare sequences? {#choose}

<a id="chord"></a>

### Chord · Overview of connections {#chord}

Show sequences and their links along curved backbones. See [backbone shape, curvature and radius](chord-geometry.md).

<a id="collinear"></a>

### Collinear · Order and local alignment {#collinear}

Arrange sequences in rows or columns to compare gene order, anchors and local regions. See [alignment and arrangement](collinear-layout.md).

<a id="radial"></a>

### Radial · Radial synteny {#radial}

Arrange sequences as spokes, keeping their own coordinates and showing syntenic relationships in the center. See [radial layout](radial-synteny.md).

<a id="sankey"></a>

### Sankey · Two-column interval relationships {#sankey}

Place sequences in query and subject columns and show observed interval relationships with ribbons. See [two-column layout](sankey.md).

<a id="circular"></a>

### Circular · Single-sequence annotation {#circular}

Arrange features, primers, restriction sites or quantitative tracks along a single circular sequence. See [single-sequence circular maps](circular.md).

<a id="composition"></a>

### Dotplot and Fusion {#composition}

[Dotplot](dotplot.md) shows paired HSP endpoints on two sequence axes. [Fusion](fusion.md) arranges independent plots and can connect them through the source intervals they display.

<a id="shared-techniques"></a>

## What do you want to adjust? {#shared-techniques}

| Goal | How to do it |
| --- | --- |
| Adjust arrow shapes and placement | [Feature shapes](feature-geometry.md) and [overlapping annotations](annotations.md) |
| Handle crowded labels | [Text fitting and repulsion](labels.md) |
| Add signals, coverage or other quantitative tracks | [Tracks](tracks.md) and [track statistics](tracks-statistics.md) |
| Set colors, legends or themes | [Layer styling](styling.md) and [scales](scale-families.md) |
| Set geometry separately by sequence or strand | [Argument formats](parameters.md) |
| Save at a target size or troubleshoot the layout | [Output](output.md) and [diagnostics](diagnostics.md) |

<a id="prepare"></a>

## Shared plotting grammar {#prepare}

| Part | Role | Interface |
| --- | --- | --- |
| Reference data | Supplies sequence IDs, lengths and units | `ggchord()` |
| Geometric layers | Draw backbones, links, features, tracks and text | `geom_*()` |
| Mapping | Connects fields to visual aesthetics | `aes()` |
| Scales | Set how values and legends appear | `scale_*()` |
| Local positions | Arrange features, labels and tracks | `position_*()` |
| Coordinates | Set the overall layout of a plot | `coord_*()` |
| Theme | Sets fonts, spacing and other non-data appearance | `theme_*()` |

<a id="terms"></a>

## Source identity and units {#terms}

Source IDs, coordinate units and full sequence lengths come from the input. A source can have multiple display instances, and display orientation can be set independently. DNA uses bp and proteins use aa. To move between them, use an [explicit coordinate projection](projection.md).

{{terms}}

<span id="check-result"></span><span id="coordinate-controls"></span>

<a id="getting-started"></a>

## Get started {#getting-started}

[Install ggchord and draw a first plot](installation.md#first-plot). To use your own tables, see the [essential fields and coordinate conventions](data.md#tables).
