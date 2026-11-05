<div class="lesson-intro">Bring four genomes into one connected view. Combine straight and curved backbones, strand-aware genes and a continuous identity scale without hiding where the data came from.</div>

## What this view shows {#overview}
{{figure:showcasechord}}

Compare the dark blue `MT118296.1` and teal `OQ646790.1` backbones first. Their radii are 0.82 and 1.18. Orange `OR222515.1` and red `MT108731.1` use 1.02 and 1.08. Radius separates backbones visually. Read sequence lengths from the bp ticks. [Compare radius on its own](chord-geometry.html#shape).

CDS are sampled separately for each available source strand, without a length threshold. Neither MT record has a negative-strand CDS. The sole negative-strand CDS on OQ646790.1 is 207 bp. An explicit Position places positive-strand arrows inside and negative-strand arrows outside. Arrow direction, colour and data retain the source strands. Positive-strand sampling favours CDS that overlap displayed link endpoints so the view demonstrates local ribbon clearance. Every arrow uses its original annotation as the label, including hypothetical protein.

## Prepare original records {#prepare}
Extract the [tutorial input bundle](../assets/tutorial-inputs.zip) into the working directory, so it contains `examples/`. Install the optional example packages with `install.packages(c("dplyr", "purrr", "tidyr"))`; they prepare tables, map lists and inspect endpoint columns.

Put the preparation script in your working directory. It contains the input, selection and display settings.

{{recipe:chord}}

The input contains four unchanged FASTA records, their GFF3 annotations and 32 original HSPs. The display keeps HSPs at least 900 bp long and at least 75% identical, with at most three per sequence pair. This reduces clutter without merging or inventing alignments.

## Compose the layers {#compose}
The endpoint-overlap selection uses dplyr ≥ 1.1.0.

{{code:showcasechord}}

| Setting | Purpose in this plot |
| --- | --- |
| `seq_orientation`, `seq_curvature`, `seq_radius` | Adjust display direction, bend and radius per sequence. Zero curvature draws a straight backbone. |
| `seq_gap = 0.06`, `rotation = 50` | Add space between backbones and turn the whole layout to position the elements. |
| `link_avoid = "smooth"`, `alpha = 0.46` | Leave smooth clearance near gene arrows and show overlaps through translucent ribbons. |
| `position = chord_position`, `width = 0.025` | Place arrows on opposite sides by source strand with offsets of ±0.030. Narrower arrows sit close to the backbones, and labels use the same position. |
| `label_layout = "radial"`, `label_wrap = 24` | Arrange gene labels radially and wrap longer names. |
| `pident` and the default stepped scale | Map percent identity to ribbon colour and show the default stepped scale. |

## Read and adapt the result {#read}
Start with the sequence colours and gene directions, then trace the larger aligned intervals. A ribbon crossing is not evidence of a rearrangement by itself. Read alignment identity using the default stepped colour scale.

Use [link filtering](ribbons.html) to change the display threshold, [labels](labels.html) to choose annotation density, or [sequence parameters](parameters.html) to change radius and curvature. For a file-by-file introduction, start with [original sequence comparisons](real-data.html).

## Save at its intended size {#save}
```r
ggsave("four-genome-chord.png", p_show_chord, width = 12, height = 10, units = "in", dpi = 144)
```
