Use a circular map to locate features along a single reference. The bundled pUC19 record contains the 2,686 bp plasmid and its annotations. Both examples start from this record.

<a id="draw"></a>

## Draw the pUC19 backbone {#draw}

Pass the sequence table to `ggchord()` and add `geom_seq()`. Set `start_angle = 90` to place the origin at twelve o’clock. Follow the ticks clockwise to increasing positions. The first and last bases meet at the top.

{{example:puc-backbone}}

<a id="adapt"></a>

## Place nine annotated intervals on the ring {#adapt}

Add the feature table explicitly. Use `position_feature_stack()` to separate overlapping intervals, and map the recorded feature colour or type when you want to distinguish annotations. The AmpR CDS spans 1,626–2,486 bp. Near the origin, several short annotations overlap, including the M13 primer-binding intervals and the multiple cloning site.

{{example:first-circular}}

Next, [label the two CDS while retaining all nine features](annotations.md), or [plot GC content along pUC19](tracks.md). The [full plasmid map](case-circular.md) marks 39 restriction-site records at 33 cut positions. The [plasmid collection](plasmid-catalog.md) includes examples from twelve other references.
