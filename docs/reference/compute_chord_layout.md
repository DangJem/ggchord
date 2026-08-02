# Compute the chord layout

Pre-computes the coordinates of all geometric elements (sequence arcs,
ribbons, gene arrows, axes, etc.) into Cartesian (x, y) coordinates and
stores them in a layout list.

## Usage

``` r
compute_chord_layout(
  seqs,
  lens,
  seq_labels,
  seq_colors,
  seqRadius,
  seq_curvature,
  orientation,
  seq_gap,
  ribbon_data = NULL,
  ribbonGap,
  ribbon_color_scheme,
  ribbon_colors,
  ribbon_alpha,
  ribbon_ctrl_point,
  gene_data = NULL,
  geneGap,
  geneWidth,
  geneLabelRadialOffset,
  geneLabelCircumOffset,
  geneLabelCircumLimit,
  geneLabelRotation,
  gene_label_show,
  gene_label_size,
  gene_label_repel = FALSE,
  gene_label_wrap = NULL,
  gene_label_max_overlaps = Inf,
  gene_label_seed = 123,
  gene_color_scheme,
  gene_colors,
  gene_order,
  seq_label_text = NULL,
  seq_label_radius = NULL,
  seq_label_rotation = NULL,
  seq_label_size = NULL,
  axisGap,
  axisMaj,
  axisMajLen,
  axisMin,
  axisMinLen,
  labelSize,
  labelOffset,
  axisLabelOrientation,
  show_axis,
  rotation,
  debug = FALSE
)
```

## Arguments

- seqs:

  Vector of sequence IDs (order already processed)

- lens:

  Named vector of sequence lengths (names = seq_id)

- seq_labels:

  Named vector of sequence labels

- seqRadius:

  Named vector of sequence radii

- seq_curvature:

  Named vector of sequence curvatures

- orientation:

  Named vector of sequence orientations (1 or -1)

- seq_gap:

  Named vector of sequence gap proportions

- ribbon_data:

  Alignment data (already validated)

- ribbonGap:

  Named vector of ribbon gaps

- gene_data:

  Gene data (already validated)

- rotation:

  Global rotation angle (degrees)

- debug:

  Whether to output debug information

## Value

A chord layout list
