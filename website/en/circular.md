Use a circular map to locate features along a single reference. The bundled pUC19 record contains the 2,686 bp plasmid and its annotations. Both examples start from this record.

<a id="draw"></a>

## Draw the pUC19 backbone {#draw}

Pass the sequence table to `ggchord()` and add `geom_seq()`. The default circular coordinates put the origin at twelve o’clock. Follow the ticks clockwise to increasing positions. The first and last bases meet at the top.

{{example:puc-backbone}}

<a id="adapt"></a>

## Place nine annotated intervals on the ring {#adapt}

Add one feature layer from the same record. The plasmid preset reads the recorded feature styles. The long green arrow lies on the arc occupied by the AmpR CDS (1,626–2,486 bp). The cluster near the origin includes the M13 primer-binding intervals and the multiple cloning site.

{{example:first-circular}}

Next, [label the two CDS while retaining all nine features](annotations.md), or [plot GC content along pUC19](tracks.md). The [full plasmid map](case-circular.md) marks 39 restriction-site records at 33 positions. The [plasmid collection](plasmid-catalog.md) includes examples from twelve other references.
