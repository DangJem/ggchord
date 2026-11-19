A dotplot puts positions on the query on the x axis and positions on the subject on the y axis. The same four HSP endpoints used for a ribbon now describe the two ends of a line segment.

<a id="draw"></a>

## Plot one observed alignment as a segment {#draw}

Select both full-length references and their first HSP. Supply the row to `geom_alignment_segment()` and use `coord_dotplot()` to assign query and subject axes. This view needs only the segment layer. The right margin leaves room for the final coordinate label.

{{example:first-dotplot}}

The line runs from (26,298, 7,121) to (51,139, 31,959) bp. Both coordinates increase, giving a rising segment. Endpoint order carries this information, so retain it when reading a table that includes reverse alignments.

<a id="adapt"></a>

## Keep the axes and add the shorter segments {#adapt}

Remove the one-row selection and keep the same axes. The two added segments lie below and to the left of the long one. At full-genome scale the 302 bp HSP near query position 20,376 bp and subject position 403 bp appears almost as a dot.

{{example:dotplot-pair}}

Next, [read a matrix of four references](case-dotplot.md). Each occupied cell contains HSPs for one ordered pair. An empty cell has no corresponding record in the supplied table. To connect different coordinate systems on one canvas, continue to [Fusion](fusion.md).
