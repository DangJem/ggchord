# Add an axis layer

Draws axes for each sequence in the chord diagram (including axis lines,
major/minor ticks, and labels). Axis parameters (spacing, tick
count/length, label size/orientation, etc.) are specified here.

## Usage

``` r
geom_axis(
  mapping = NULL,
  data = NULL,
  show_axis = NULL,
  axis_gap = NULL,
  axis_tick_major_number = NULL,
  axis_tick_major_length = NULL,
  axis_tick_minor_number = NULL,
  axis_tick_minor_length = NULL,
  axis_label_size = NULL,
  axis_label_offset = NULL,
  axis_label_orientation = NULL,
  show_legend = FALSE,
  ...
)
```

## Arguments

- mapping:

  Default NULL (uses pre-computed data)

- data:

  Default NULL (retrieved automatically from the layout)

- show_axis:

  Logical. Whether to show the axis, default TRUE

- axis_gap:

  Optional numeric/vector. Spacing between sequence and axis, default
  0.04

- axis_tick_major_number:

  Optional integer/vector. Number of major ticks, default 5

- axis_tick_major_length:

  Optional numeric/vector. Major tick length ratio, default 0.02

- axis_tick_minor_number:

  Optional integer/vector. Number of minor ticks, default 4

- axis_tick_minor_length:

  Optional numeric/vector. Minor tick length ratio, default 0.01

- axis_label_size:

  Optional numeric/vector. Tick label font size, default 3

- axis_label_offset:

  Optional numeric/vector. Label offset ratio, default 1.5

- axis_label_orientation:

  Optional character/numeric/vector. Label orientation, default
  "horizontal"

- show_legend:

  Whether to show the legend, default FALSE (axes do not participate in
  legends)

- ...:

  Additional arguments passed to geom_path/geom_segment/geom_text

## Value

A list of ggplot2 layers
