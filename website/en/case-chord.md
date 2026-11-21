Four assembled sequences differ in length and share several aligned regions. This map separates the backbones with a mixture of straight and curved segments so that the CDS and ribbons remain easy to locate.

<a id="map"></a>

## Draw the full map {#map}

The figure uses a 12 × 10 in canvas at 144 dpi. The script returns the plot so you can pass it to `ggplot2::ggsave()` at the same size.

{{example:case-chord}}

<a id="reading"></a>

## Read the figure {#reading}

Start at a sequence label. The backbone scale gives position in bp. A ribbon encloses an observed alignment interval at each end, and its fill encodes percent identity. Arrows on either side of the backbone mark CDS with a recorded biological strand.

<a id="inputs"></a>

## Source records and display choices {#inputs}

`flagship_comparison_example$chord` holds the four lengths, 12 HSPs and 20 displayed CDS. The selection kept alignments at least 900 bp long with at least 75% identity, then retained up to three of the longest per sequence pair. Named settings in `layout_sequences()` control display orientation, radius and curvature. Each setting matches a reference accession.

[Draw the first ribbon](chord.md) · [See the same input on radial spokes](case-radial.md)
