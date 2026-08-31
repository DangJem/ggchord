# Genomic sequence position scale

Controls major/minor genomic position breaks and labels independently
for every sequence. The scale is trained against each sequence's \`\[0,
length\]\` range during chord layout.

## Usage

``` r
scale_seq_position_continuous(
  name = ggplot2::waiver(),
  breaks = ggplot2::waiver(),
  minor_breaks = ggplot2::waiver(),
  labels = ggplot2::waiver(),
  limits = NULL,
  expand = ggplot2::waiver(),
  oob = scales::censor,
  transform = "identity",
  ...
)
```

## Arguments

- name:

  Scale name; position guides are disabled by default.

- breaks, minor_breaks, labels, limits, expand, oob, transform:

  Standard continuous-scale controls.

- ...:

  Additional arguments passed to \[ggplot2::continuous_scale()\].

## Value

A ggplot2 continuous scale for the \`seq_position\` role.
