The radial map places four full-length sequences on separate spokes. It uses the same 12 HSPs and 20 CDS as the chord case, so the two geometries can be compared directly.

<a id="map"></a>

## Draw the full map {#map}

The figure uses a 12 × 10.5 in canvas at 144 dpi. The script returns the plot so you can pass it to `ggsave()` at the same size.

{{example:case-radial}}

<a id="reading"></a>

## Read the figure {#reading}

A coloured sequence backbone carries its own bp scale. The arrow and label beside it locate a CDS. For ribbons between non-neighbouring spokes, the outward bend gives the connection room to pass around intermediate spokes.

<a id="inputs"></a>

## Source records and display choices {#inputs}

Use the stored orientation vector in `layout_sequences()`, which also controls the size of the central opening. Set `start_angle = 45` in `coord_radial_synteny()`. Ribbon curvature and the outward bend are controlled on the ribbon layer. The recorded endpoints remain in source coordinates.

[First radial comparison](radial-synteny.md) · [Chord view of these sequences](case-chord.md)
