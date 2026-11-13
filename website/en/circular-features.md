{{code:courseblockak}}

{{prepare:comparefinishcircle}}

## Generic shapes and the plasmid preset {#feature-shapes}

On the same six source annotations, change only arrow to block. A block omits the directional arrow without changing strand. The third figure shows the combined geom_feature_plasmid preset: width=0.058, arrow_head_length=0.068, arrow_head_width=1.42 and plasmid Position. It changes several defaults and is not a single-parameter contrast. See [shape choices](annotations.html#shapes) and [arrow/short-feature geometry](feature-geometry.html).

{{code:comparefeaturearrow}}
{{figure:comparefeaturearrow}}
{{code:comparefeatureblock}}
{{figure:comparefeatureblock}}
{{code:comparefeaturepreset}}
{{figure:comparefeaturepreset}}
## Share placement between overlapping features and labels {#stacking}

Compare identity versus a shared feature stack on the same six annotations. Pass the Position to both geometry and labels rather than moving text alone. `position_feature_stack(spacing = 0.1, side = "strand", base_position = NULL)` applies an optional base Position before allocating overlap lanes. `position_plasmid()` is a local preset; `position_strand()` applies signed offsets by biological strand. Inside/outside is not a substitute for strand evidence.

{{code:comparefinishstack}}
{{code:comparepositionidentity}}
{{figure:comparepositionidentity}}
{{code:comparepositionstack}}
{{figure:comparepositionstack}}
## Exercises, replacement input and next steps {#adapt}

Change only double-line spacing to 0.04, then change the centre length_style size. Keep complete source lengths; do not infer them from maximum feature ends. For FASTA/GenBank/GFF3, check IDs and coordinates in [formats](formats.html). Continue with [circular labels](circular-labels.html), [primers](primers.html), [restriction sites](restriction-sites.html) and [tracks](tracks.html). Preview at the figure dimensions and use the same width/height in ggsave.
