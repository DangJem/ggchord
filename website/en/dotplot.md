A dotplot maps query coordinates to x and subject coordinates to y. Each alignment segment connects two recorded coordinate pairs. It shows which regions align, their spans and their relative direction.

<a id="draw"></a>

## Draw one observed segment {#draw}

This segment starts at (26,298, 7,121) bp and ends at (51,139, 31,959) bp. Both source coordinates increase between the endpoints, so the segment has positive slope on these axes.

{{example:first-dotplot}}

<a id="adapt"></a>

## Read a matrix of sequence pairs {#adapt}

Supply full sequence lengths and retain the order of both endpoint pairs. `x_seqs` and `y_seqs` select and order the two axes independently. In a multi-sequence matrix, an empty cell contains no record from the supplied table. Per-base gaps and mismatches require alignment columns or a finer alignment path.

[Four-sequence alignment matrix](case-dotplot.md) · [Alignment segment layer](../../man/geom_alignment_segment.Rd) · [Dotplot axes](../../man/coord_dotplot.Rd)
