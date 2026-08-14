# Add a sequence arc layer

Draws arcs (or straight lines, depending on the curvature setting)
representing sequences in the chord diagram. Sequence layout parameters
(order, orientation, radius, curvature, colors, etc.) are specified
here. Sequences can be visually grouped with an extra inter-group gap
and optional group labels.

## Usage

``` r
geom_seq(
  mapping = NULL,
  data = NULL,
  seq_order = NULL,
  seq_labels = NULL,
  seq_orientation = NULL,
  seq_gap = NULL,
  seq_radius = NULL,
  seq_curvature = NULL,
  seq_colors = NULL,
  seq_group = NULL,
  seq_group_gap = 0.08,
  seq_group_labels = TRUE,
  seq_group_label_radius = 1.35,
  seq_group_colors = NULL,
  linewidth = 1.2,
  show_legend = TRUE,
  legend_position = "right",
  ...
)
```

## Arguments

- mapping:

  Default NULL (uses pre-computed data)

- data:

  Default NULL (retrieved automatically from the layout)

- seq_order:

  Optional character vector. Specifies the drawing order of sequences

- seq_labels:

  Optional character vector or named vector. Sequence labels

- seq_orientation:

  Optional numeric (1 or -1). Sequence orientation, default 1

- seq_gap:

  Optional numeric. Gap proportion between sequences, default 0.03

- seq_radius:

  Optional numeric (\> 0). Sequence arc radius, default 1.0

- seq_curvature:

  Optional numeric. Arc curvature (0=straight, 1=standard arc, \>1=more
  curved), default 1.0

- seq_colors:

  Optional color vector or named vector. Sequence colors

- seq_group:

  Optional group specification. NULL disables grouping unless `seq_data`
  contains a `seq_group` column; a single string naming such a column
  uses that column; otherwise a single value, named vector, unnamed
  vector matching the sequences, or a list (see
  [`process_sequence_param()`](https://dangjem.github.io/ggchord/reference/process_sequence_param.md))
  is accepted.

- seq_group_gap:

  Numeric, default 0.08. Extra gap proportion inserted between
  consecutive groups (in addition to `seq_gap`).

- seq_group_labels:

  Logical or character, default TRUE. When TRUE the group names are
  drawn at the angular midpoint of each group; a named character vector
  can be used to override the group label text.

- seq_group_label_radius:

  Numeric, default 1.35. Radial position of the group labels as a
  multiplier of the group's outermost sequence radius (same convention
  as `seq_label_radius`: `1` on the arc, `> 1` outside).

- seq_group_colors:

  Optional named color vector by group name (or an unnamed vector
  recycled positionally). Colours the group labels.

- linewidth:

  Arc line width, default 1.2

- show_legend:

  Whether to show the legend for this layer, default TRUE

- legend_position:

  Position of this layer's legend (the Seq ID legend): one of "left",
  "right", "top", "bottom" or "inside", default "right". Pass NULL to
  let the legend follow `theme(legend.position = ...)` together with the
  other legends. Can also be set with
  `theme(legend.position.seq = ...)`.

- ...:

  Additional arguments passed to `geom_path()`

## Value

A list of ggplot2 layers

## Examples

``` r
library(ggchord)
data(seq_data_example)
p <- ggchord(seq_data_example) + geom_seq()
p
```
