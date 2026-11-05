<div class="lesson-intro">The complete Module II result: four real genomes, strand-separated features and outward-routed cross-row links.</div>

## Complete result {#overview}

{{figure:showcasecollinear}}

Compare four backbones on a shared base-pair scale with 16 selected CDS placed on either side of the rows. Read gene positions along one row first, then follow a band to its matching interval. A solid outline marks adjacent rows, while a dashed outline marks rows with another displayed row between them. Start with the [four-backbone plot](collinear-start.html) and [alignment and CDS example](collinear-build.html) if you are new to the workflow.

## Prepare original records {#prepare}
Extract the [input bundle](../assets/tutorial-inputs.zip) and install dplyr, purrr and tidyr first:

Put the preparation script in your working directory. It contains the input, selection and display settings.

{{recipe:collinear}}

Retain four full FASTA records, 266 CDS and 32 original HSPs. Adjacent rows use the explicitly selected observations; each non-neighbour pair retains only its longest HSP to reduce overlap.

## Reserve space for connections and labels {#gc}
To compare sequence composition, see [tracks](tracks.html) for a GC track example.

## Compose the layers {#compose}
{{code:showcasecollinear}}

Dashed outlines identify non-adjacent displayed rows, with sufficient linewidth to read. Right-side controls separate long routes; smooth local clearance retreats near arrows. Identity has the same left-side colourbar as the chord view, and the bottom strand keys use arrows. The ruler's horizontal offset is 0; its length still uses the shared bp conversion.

## Adjustment and limits {#read}
Rebuild controls when changing seq_gap, order, orientation or link selection. Local arrow clearance is not an automatic route around intervening whole sequences; ribbons may still cross legitimately. Use [focus/segments](views.html) for local questions. Explicit circular radius/curvature cannot currently be passed to linear coordinates; see [conversion boundaries](collinear-data.html#contract).

## Save at the intended size {#save}
```r
ggsave("four-genome-collinear.png", p_show_collinear, width = 12, height = 8, units = "in", dpi = 144)
```
