<div class="task-intro"><p>Reduce link clutter while retaining provenance: filtering, deduplication, merging and bundling have different meanings.</p></div>

## Prepare this page {#prepare}
In a fresh R session, run the code below, then follow this page in order. Colours and layout are optional display choices; inputs select only the fields needed. For file-based examples, first extract the [input bundle](../assets/tutorial-inputs.zip) so the working directory contains `examples/`.

{{prepare:setup,docstyle,realimport,realprepare}}

## Begin with the original link table {#input}
Use the [real-record workflow](real-data.html#selection) to create `seq_real`,
`link_real` and `real_filtered`. Keep its original tables and filter report.

{{code:ribbons}}
{{results:ribbons}}

## Understand each transformation {#transformations}
Filtering keeps observations meeting declared thresholds. Reciprocal-overlap
deduplication selects a representative existing observation; it does not sum
identities. Merging combines nearby compatible intervals into a display envelope.
Bundling groups nearby query/subject locations into an aggregate. Inspect each
`$report` and source-row membership, especially when extra metadata conflict.
The merged or bundled envelope is not a newly computed nucleotide alignment.
Estimated lengths retain their flags; unobserved identity is never invented.

## Density is an explicit statistic {#density}
{{figure:density}}

`stat_ribbon_density()` maps a computed density, and `stat_ribbon_bundle()` performs
its explicit aggregation during build. The ordinary ribbon geom still draws one
ribbon per row. Interpret count, reported length and identity weighting in their
own units; more opaque ribbons need not mean stronger independently validated
homology. A point anchor without direction evidence cannot acquire an aggregate
orientation merely because its endpoints are equal.

## Choose interval ribbons or point links {#points}
After this page’s preparation code:

{{code:pointlinks}}
{{figure:pointlinks}}

Point positions are integer midpoints of real bundled HSPs; identity remains the reported `pident`. Lines use `qpos/spos`,
not interval ends, and an arrow points from query to subject. The line default
avoidance is `none`; ribbon avoidance is `uniform`. Global chord rotation belongs
to the coordinate system. Pair selection and adjacent-row defaults in linear
plots are covered in [three modes](modes.html#collinear).

The point-link colour scale supplies `available_aes = "link_colour"` to its
colourbar explicitly. Ordinary ggplot2 colourbars and the wrapper default do not
accept every extended role; choose the matching role when making a custom guide.

## Check your result {#check-result}
Inspect each report and member rows before reading the plot. Aggregated bands are not new alignments; point links summarize midpoints of existing HSPs.

<div class="next-step">Continue: <a href="views.html">Focus, segments and homology →</a></div>
