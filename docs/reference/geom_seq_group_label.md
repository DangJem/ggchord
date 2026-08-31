# Add sequence-group labels

Draws labels for contiguous sequence groups defined by \`seq_group\` in
\`geom_seq()\`. This dedicated layer lets group-label typography be
controlled independently through \`ggchord.group.label\`.

## Usage

``` r
geom_seq_group_label(
  mapping = NULL,
  data = NULL,
  labels = TRUE,
  radius = 1.35,
  size = NULL,
  show_legend = FALSE,
  ...
)
```

## Arguments

- mapping, data:

  Optional layer mapping and sequence data.

- labels:

  Logical or character labels. \`TRUE\` uses group names; named values
  replace selected group names.

- radius:

  Radial position relative to the outermost sequence radius.

- size:

  Optional text size in millimetres. The theme element is used when
  \`NULL\`.

- show_legend:

  Whether to show a group-colour legend.

- ...:

  Additional fixed arguments passed to \`geom_text()\`.

## Value

A list containing one ggplot2 layer.
