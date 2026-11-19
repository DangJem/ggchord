{{code:courseblockak}}


<a id="feature-shapes"></a>

## Generic shapes and the plasmid preset {#feature-shapes}

On the same six source annotations, change only arrow to block. A block omits the directional arrow without changing strand. The third figure shows the combined geom_feature(preset = 'plasmid') preset: width=0.058, arrow_head_length=0.068, arrow_head_width=1.42 and plasmid Position. It changes several defaults and is not a single-parameter contrast. See [shape choices](annotations.md#shapes) and [arrow/short-feature geometry](feature-geometry.md).

{{code:comparefeaturearrow}}
{{figure:comparefeaturearrow}}
{{code:comparefeatureblock}}
{{figure:comparefeatureblock}}
{{code:comparefeaturepreset}}
{{figure:comparefeaturepreset}}

<a id="stacking"></a>

## Share placement between overlapping features and labels {#stacking}

Compare identity and shared stacking on the same six annotations. Pass the same Position to geometry and text so their lanes agree. `position_feature_stack(spacing = 0.1, base_position = position_strand(), direction = "auto")` places the starting lane before allocating overlaps. Use `position_offset()` as a shared starting point for both strands, with the growth direction set independently. Both layouts preserve source strand.

{{code:comparefinishstack}}
{{code:comparepositionidentity}}
{{figure:comparepositionidentity}}
{{code:comparepositionstack}}
{{figure:comparepositionstack}}

<a id="adapt"></a>

## Replace input and continue {#adapt}

Change only double-line spacing to 0.04, then change the centre length_style size. Keep complete source lengths; do not infer them from maximum feature ends. For FASTA/GenBank/GFF3, check IDs and coordinates in [formats](formats.md). Continue with [circular labels](circular-labels.md), [primers](primers.md), [restriction sites](restriction-sites.md) and [tracks](tracks.md). Preview at the figure dimensions and use the same width/height in ggsave.
