A circular map places positions from one reference around a ring. Start with its full length and an annotation table. Feature coordinates share the backbone’s units, so the positions at the beginning and end meet at the seam.

<a id="draw"></a>

## Draw the nine pUC19 features {#draw}

The installed pUC19 record provides the 2,686 bp length and nine existing annotations. `geom_feature()` reads the feature table directly. The arrows follow the annotated strand where a direction is known.

{{example:first-circular}}

<a id="adapt"></a>

## Add labels and sites {#adapt}

Attach the same position adjustment to `geom_feature()` and `geom_feature_label()` to keep labels with their lanes. Use separate layers for primer intervals and restriction sites. The full pUC19 record includes 39 site records at 33 distinct positions. A shared position can carry more than one enzyme name.

[Fully annotated pUC19 map](case-circular.md) · [Thirteen plasmid maps](plasmid-catalog.md) · [Circular coordinates](../../man/coord_circular.Rd)
