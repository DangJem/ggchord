Sequence names identify entire backbones. Gene names identify individual intervals. Their label layers can be adjusted separately.

## Move sequence names inward {#sequence-labels}

{{code:frontseqlabels}}

{{figure:frontseqlabels}}

`seq_label_radius = 0.7` places each label at 70% of its backbone's radius. Compare with the [default plot](chord-start.html#draw). The labels sit farther from the ticks. Text follows the arc by default. `seq_label_orientation = "horizontal"` makes all four names horizontal. Change one setting at a time to separate their effects.

## Label one CDS {#gene-labels}

{{code:frontgenes,frontgenelabels}}

{{figure:frontgenelabels}}

Find the red integrase arrow on the green `OQ646790.1` backbone. The new label comes from `anno`. Here, `label_radial_offset = 0.10` puts a little more radial distance between the text and arrow. When several annotations are close together, use `geom_gene_label_repel()` to spread their text. See [dense labels](labels.html#layout).

Continue to the [four-sequence case study](chord.html) to combine these adjustments.
