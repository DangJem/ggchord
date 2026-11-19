To compare positions along two genomes, lay their backbones on parallel rows. Start with the two phages and the single local alignment (HSP) from the [chord diagram](chord.md).

<a id="draw"></a>

## Locate one alignment on two rows {#draw}

Select the first HSP from MT108731.1 to MT118296.1. The ribbon layer reads its four endpoints. Changing the last line to `coord_collinear()` gives each reference its own horizontal axis in bp.

{{example:views-collinear}}

The upper interval runs from 26,298 to 51,139 bp. The lower one runs from 7,121 to 31,959 bp and nearly reaches the end of the 32,090 bp subject. The shorter lower backbone shows the difference in genome lengths.

<a id="adapt"></a>

## Keep every HSP for this ordered pair {#adapt}

Keep the two accession filters and drop the one-row selection. Two shorter ribbons now lie to the left of the long ribbon. The narrowest links 20,376–20,677 bp on the query to 403–704 bp on the subject.

{{example:first-collinear}}

Next, [read the first alignment along radial spokes](radial-synteny.md). To reverse a row or repeat a reference, try [sequence layout](multi-sequence-layout.md). The [four-genome case](case-collinear.md) adds CDS annotations and nonadjacent comparisons.
