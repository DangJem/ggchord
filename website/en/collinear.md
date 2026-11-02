<div class="lesson-intro">The complete Module II result: four real genomes, strand-separated features and outward-routed cross-row links.</div>

## Complete result {#overview}
{{figure:showcasecollinear}}

Collinearity extends the chord coordinate grammar; reusing data and layers is the goal. Read [inputs](collinear-data.html), [coordinate replacement](collinear-build.html) and [routes/ruler](collinear-routing.html), then run the complete example. The identity scale and feature candidates are shared with the chord example; this linear composition omits two crowded short candidates, displaying 16 CDS with labels, with positive and negative strands above/below backbones.

## Prepare original records {#prepare}
Extract the [input bundle](../assets/tutorial-inputs.zip) and install dplyr, purrr and tidyr first:

{{prepare:sequencepalette,docstyle}}
{{code:showcasedata}}

Retain four full FASTA records, 266 CDS and 32 original HSPs. Adjacent rows use the explicitly selected observations; each non-neighbour pair retains only its longest HSP to reduce overlap.

## Reserve space for connections and labels {#gc}
The flagship no longer includes GC tracks; composition analysis remains in [tracks](tracks.html). This keeps two different reading tasks from competing in the entry figure.

## Compose the layers {#compose}
{{code:showcasecollinear}}

Dashed outlines identify non-adjacent displayed rows, with sufficient linewidth to read. Right-side controls separate long routes; smooth local clearance retreats near arrows. Identity has the same left-side colourbar as the chord view, and the bottom strand keys use arrows. The ruler's horizontal offset is 0; its length still uses the shared bp conversion.

## Adjustment and limits {#read}
Rebuild controls when changing seq_gap, order, orientation or link selection. Local arrow clearance is not an automatic route around intervening whole sequences; ribbons may still cross legitimately. Use [focus/segments](views.html) for local questions. Explicit circular radius/curvature cannot currently be passed to linear coordinates; see [conversion boundaries](collinear-data.html#contract).

## Save at the intended size {#save}
```r
ggsave("four-genome-collinear.png", p_show_collinear, width = 12, height = 8, units = "in", dpi = 144)
```
