Four rows show where aligned regions fall along each full-length sequence. Plotting CDS beside the backbones helps relate these regions to the order and orientation of annotated genes.

<a id="map"></a>

## Draw the full map {#map}

The figure uses a 12 × 8 in canvas at 144 dpi. The script returns the plot so you can pass it to `ggplot2::ggsave()` at the same size.

{{example:case-collinear}}

<a id="reading"></a>

## Read the figure {#reading}

Seven selected HSPs join adjacent rows. A separate layer draws the three longer-range links across rows. Follow the two ends of a ribbon to compare source coordinates, then use the nearest CDS arrow and label to identify a feature.

<a id="inputs"></a>

## Source records and display choices {#inputs}

The collinear subset contains 16 CDS. The two link tables preserve the same six endpoint fields, so they can use independent styling and routing. If the row order or spacing changes, review the paths of the cross-row ribbons and check the same two source endpoints.

[Compare two rows](collinear.md) · [Read the same intervals in two columns](case-sankey.md)
