Complete each small step before the flagship example. All four views reuse the same records, making the effect of each added layer visible.

## Prepare this page {#prepare}
Extract the [input bundle](../assets/tutorial-inputs.zip) first.
{{prepare:sequencepalette,docstyle,showcasedata}}

## Step 1: sequence backbones {#backbones}
{{code:learnseq}}
{{figure:learnseq}}
`geom_seq()` draws backbones, `geom_seq_label()` identifies sequences, and `coord_chord()` selects the global view. No links or genes are required yet.

## Step 2: observed comparisons {#links}
{{code:learnlinks}}
{{figure:learnlinks}}
Identity belongs in aes and a scale. Explicitly retaining non-neighbour links keeps the same link data when replacing coordinates later.

## Step 3: selected genes {#genes}
{{code:learngenes}}
{{figure:learngenes}}
`position_strand()` separates positive and negative strands. Display selected real annotations while retaining every source CDS separately.

## Step 4: label every displayed feature {#labels}
{{code:learnlabels}}
{{figure:learnlabels}}
The arrow and repel-label layers receive the same feature table without another label filter. Renaming, deleting or reordering records can change meaning; retain clear intermediate objects.

Continue to [1.4 Direction, radius, curvature and clearance](chord-tuning.html).
