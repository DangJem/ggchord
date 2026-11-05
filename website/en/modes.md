## What makes up a plot? {#prepare}

| Concept | Meaning in the plot | In code |
| --- | --- | --- |
| Sequence backbone | Represents one sequence. Intervals are positioned along it. | `geom_seq()` |
| Feature | An interval on a sequence, such as a CDS or replication origin. A gene with a known strand can appear as an arrow. | `geom_feature()`, `geom_gene()` |
| Link | Joins corresponding intervals with a ribbon or joins two positions with a line. | `geom_link_ribbon()`, `geom_link_line()` |
| Layer | Adds one kind of element, such as a backbone or label, to the plot. | Add `geom_*()` with `+` |
| Coordinate system | Arranges backbones as arcs, rows or a single ring. | `coord_*()` |
| Mapping and scale | Connect table values to properties such as colour and specify how they are displayed. | `aes()` and `scale_*()` |

Biological strand comes from the feature table. Display orientation determines which way coordinates increase on the page. Alignment direction comes from endpoint order or alignment records. Identify which direction you are reading when looking at an arrow, backbone or ribbon.

## Choose a layout {#choose}

### Chord plot {#chord}

{{figure:frontbackbone}}

`coord_chord()` arranges multiple sequences in one view so you can see their overall relationships. [Start with sequence backbones](chord-start.html).

### Collinear plot {#collinear}

{{figure:frontlinearbase}}

`coord_collinear()` lays backbones out in rows to make position and order easier to compare. [Open the collinear chapter](collinear-start.html). Collinear ribbons default to links between adjacent display rows.

### Single-sequence circular plot {#circular}

{{figure:firstcircle}}

`coord_circular()` draws a single sequence as a ring, suitable for annotating plasmids and other circular DNA. [Open the circular chapter](plasmid-start.html).

## From data to a figure {#composition}

Start with full sequence IDs, lengths and units, then distinguish annotations from relationships between sequences. The same input can use different coordinate layouts. Source positions and biological strands remain the basis for interpretation.

Colour may identify sequences, feature classes or alignment identity. Legends should explain that meaning. Biological strand, display orientation and relative alignment direction have distinct sources and should be checked separately.

## Terms and units {#terms}

{{terms}}

<span id="check-result"></span>Choose a plotting module and start with default backbones. Explore parameters, labels and export settings once you have a figure to compare.

## Fit, limits and clipping {#coordinate-controls}

These controls act on the complete view after sequence geometry. They do not change source coordinates or infer missing lengths. See the [parameter atlas](parameter-atlas.html) for each coordinate signature.

| Control | Default | Practical effect |
| --- | --- | --- |
| `fit` | `"labels"` | `labels` includes measured labels; `geometry` fits geometry only; `manual` requires xlim and ylim. Recheck after changing label or export size |
| `xlim` / `ylim` | `NULL` | Explicit user limits take priority over automatic fitting; they can omit content |
| `ratio` | `1` | Fixed y/x aspect; changing it affects visual geometry, not biological distance |
| `clip` | `"off"` | `"on"` clips drawing to the panel; it does not remove source records |
| `expand` | `FALSE` | Additional ggplot2 expansion; automatic fitting already adds a small margin |

A fitted panel cannot guarantee that every label remains readable. Inspect the final target size and diagnostics; avoid solving a source-coordinate mistake by widening limits.
