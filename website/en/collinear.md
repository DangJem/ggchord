Collinear coordinates lay sequences on separate rows. Use them to compare the order and orientation of aligned intervals, or to see where an annotated feature falls relative to an alignment.

<a id="draw"></a>

## Compare two full-length sequences {#draw}

The sequence table provides full lengths in bp. The filter keeps every input alignment in the displayed query-to-subject direction. Three ribbons remain. Read their two widths against the two backbones.

{{example:first-collinear}}

<a id="adapt"></a>

## Choose row order and orientation {#adapt}

With several rows, the default ribbon layer shows adjacent pairs. Use `adjacent_only = FALSE` to include selected non-adjacent pairs. Set the order with `layout_sequences(order = ...)`. A named orientation vector can flip a displayed row while retaining its source coordinates and strand.

[Labelled four-genome comparison](case-collinear.md) · [Linear coordinates and anchoring](../../man/coord_collinear.Rd) · [Sequence ordering](../../man/layout_sequences.Rd)
