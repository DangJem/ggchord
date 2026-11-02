<div class="lesson-intro">Bring four genomes into one connected view. Combine straight and curved backbones, strand-aware genes and a continuous identity scale without hiding where the data came from.</div>

## What this view shows {#overview}
{{figure:showcasechord}}

This is the end of Module I. Begin with [tables](chord-data.html), [import](chord-import.html), [layer construction](chord-build.html) and [geometry/clearance](chord-tuning.html). 18 selected real CDS each receive repel labels; the complete 266-row annotation source remains separate. The four backbones combine zero, positive and negative curvature, varying radii and optimized display directions.

## Prepare original records {#prepare}
Extract the [tutorial input bundle](../assets/tutorial-inputs.zip) into the working directory, so it contains `examples/`. Install the optional example packages with `install.packages(c("dplyr", "purrr", "tidyr"))`; they prepare tables, map lists and inspect endpoint columns.

{{prepare:sequencepalette,docstyle}}
{{code:showcasedata}}

The input contains four unchanged FASTA records, their GFF3 annotations and 32 original HSPs. The display keeps HSPs at least 900 bp long and at least 75% identical, with at most three per sequence pair. This reduces clutter without merging or inventing alignments.

## Compose the layers {#compose}
{{code:showcasechord}}

`geom_seq(seq_curvature = ...)` applies named values to stable sequence identifiers. `0` is a straight backbone. Curvature changes the display geometry; it does not alter sequence lengths, biological strands or alignment endpoints.

## Read and adapt the result {#read}
Start with the sequence colours and gene directions, then trace the larger aligned intervals. A ribbon crossing is not evidence of a rearrangement by itself. The colour bar uses the same identity scale as the [collinear view](collinear.html).

Use [link filtering](ribbons.html) to change the display threshold, [labels](labels.html) to choose annotation density, or [sequence parameters](parameters.html) to change radius and curvature. For a file-by-file introduction, start with [original sequence comparisons](real-data.html).

## Save at its intended size {#save}
```r
ggsave("four-genome-chord.png", p_show_chord, width = 12, height = 10, units = "in", dpi = 144)
```
