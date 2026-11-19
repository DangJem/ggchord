Sankey coordinates display the references in two columns. A ribbon starts at a query interval on the left and ends at its subject interval on the right. Its width at each end follows the recorded span.

<a id="draw"></a>

## Follow a query across columns {#draw}

The first three input rows link MT108731.1 to MT118296.1. Locate the three intervals on the left copy of the first accession, then follow each band to the right copy of the second. Each reference appears in both columns.

{{example:first-sankey}}

<a id="adapt"></a>

## Interpret width and direction {#adapt}

Width is measured in the coordinate units of the endpoints. An alignment can span different numbers of positions at its two ends. Reciprocal query-to-subject rows count as separate observations when supplied. Adjust `column_gap` to separate the columns.

[Four-genome two-column map](case-sankey.md) · [Sankey coordinates](../../man/coord_sankey.Rd)
