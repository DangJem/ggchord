<a id="anchors"></a>

## Draw all features, label selected records {#anchors}

The pUC19 example contains nine reference features. A stable ID links each label to its source feature. Allocate lanes using the full table, then select the two CDS records for labels. The entity layer has the ID `features`, which the label layer names in `anchor_layer`.

{{example:feature-subset-labels}}

<a id="placement"></a>

## Choose what controls the geometry {#placement}

Use `geom_feature()` for intervals with mixed or unknown strand. The gene preset draws directed coding features on nucleotide references with known strand. Supply the source strand when a teaching table also contains a recoded display strand. For same-side lanes, give `position_feature_stack()` a `base_position = position_offset(...)` and a lane direction. The feature key and original interval let subset labels follow those lanes after a reorder.

For a label drawn at a chosen angle, use fixed placement. For automatic placement, use `placement = "repel"`. Choose `overflow` for text that exceeds the available space. Inspect labels at the [intended export dimensions](output.md).
