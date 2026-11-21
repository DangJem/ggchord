Four source sequences appear in both columns. This arrangement lets the eye follow the query end on the left and the subject end on the right even when the two accessions are far apart in the row order.

<a id="map"></a>

## Draw the full map {#map}

The figure uses a 10 × 16 in canvas at 144 dpi. The script returns the plot so you can pass it to `ggplot2::ggsave()` at the same size.

{{example:case-sankey}}

<a id="reading"></a>

## Read the figure {#reading}

The input includes ten selected HSPs and 16 CDS. Percent identity controls ribbon fill. Arrows mark CDS strand and the sequence colours identify repeated copies of each accession. A second copy of a sequence adds a display location for its features.

<a id="inputs"></a>

## Source records and display choices {#inputs}

The adjacent and cross-row tables come from the collinear display subset. `layout_sequences()` sets the column gap in millimetres. Both columns share the sequence and strand legends. For a narrower export, first check labels on the outer sides.

[Read a two-column diagram](sankey.md) · [Compare the same rows in a linear layout](case-collinear.md)
