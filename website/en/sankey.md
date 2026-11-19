Use a two-column Sankey view when you want to distinguish the query and subject ends of an alignment. Each reference appears in both columns, and each input row connects a left query interval to a right subject interval.

<a id="draw"></a>

## Trace one band across the two columns {#draw}

Keep the two phages and the first HSP. Read the band from MT108731.1 on the upper left to MT118296.1 on the lower right. The two other copies provide coordinate context and carry no endpoint for this record.

{{example:views-sankey}}

<a id="adapt"></a>

## Compare the widths of three HSPs {#adapt}

Include all three rows in the same query-to-subject direction. Each band takes its width from the span on that particular reference. The broad band covers roughly 25 kb on each side. The shortest covers just 302 bp on each side.

{{example:first-sankey}}

A reciprocal row in the input would give a separate band in the opposite sequence pairing. For a denser layout, `column_gap` sets the distance between columns. Next, [compare all four layouts](multi-sequence-views.md) or [plot the HSP as an x–y segment](dotplot.md). The [four-genome Sankey case](case-sankey.md) adds annotations to the repeated references.
