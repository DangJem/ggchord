This is the full-size pUC19 map linked from the circular preview. It overlays the saved source annotations and restriction positions on a 2,686 bp double backbone.

<a id="map"></a>

## Draw the full map {#map}

The figure uses a 10 × 8 in canvas at 110 dpi. The script returns the plot so you can pass it to `ggsave()` at the same size.

{{example:case-circular}}

<a id="reading"></a>

## Read the figure {#reading}

The nine feature intervals retain the reference coordinates and strand. Follow the backbone ticks to read the features, then look outside the circle for the restriction-site names. The 47 site records occupy 32 positions. Enzymes that share a position have a combined label. This record contains no independent reference primers.

<a id="inputs"></a>

## Source records and display choices {#inputs}

`flagship_circular_example$references[["pUC19"]]` contains the complete length and the three stored plotting tables. The feature layer sets the stacking position, and the label layer follows it by ID. The code adds primer geometry only for a reference that has primer rows.

[Draw only the nine features](circular.md) · [Compare all thirteen reference maps](plasmid-catalog.md)
