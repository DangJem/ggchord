On a dense circular map, the set of features you want to draw is often larger than the set you want to label. Keep the whole feature table for placement, then select labels by a recorded field.

<a id="placement"></a>

## Locate the features before choosing labels {#placement}

Draw the nine pUC19 annotations on their 2,686 bp backbone. Each feature row carries an accession, start and end, and a strand where known. In the cluster near the origin, several short intervals occupy overlapping spans.

{{example:first-circular}}

<a id="anchors"></a>

## Label two CDS while placing all nine features {#anchors}

Keep the source match ID as `feature_id`. The layer named `features` places overlapping intervals in lanes with `position_feature_stack()`. The label layer selects the two CDS and names the original layer in `anchor_layer` so each label follows its feature, even after the selection. `layout_label_repel()` arranges the text around those anchors. Set `leader` and `segment_colour` on the label layer to control the connecting lines.

{{example:feature-subset-labels}}

Look for lacZα close to the origin and AmpR across the longer interval on the opposite side of the map. The seven unlabelled records still appear. Next, [map annotation types to a palette](styling.md), or [see how the full pUC19 map labels sites](case-circular.md).
