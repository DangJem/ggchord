<div class="task-intro"><p>Reduce link clutter while retaining provenance: filtering, deduplication, merging and bundling have different meanings.</p></div>

<span id="prepare"></span>

## Begin with the original link table {#input}

{{code:ribboninputs,ribbons}}

{{results:ribbons}}

The first three counts follow one pipeline. Bundling starts from the original table and has its own count.

## Understand each transformation {#transformations}
Filtering keeps observations meeting declared thresholds. Reciprocal-overlap
deduplication selects a representative existing observation; it does not sum
identities. Merging combines nearby compatible intervals into a display envelope.
Bundling groups nearby query/subject locations into an aggregate. Inspect each
`$report` and source-row membership, especially when extra metadata conflict.
The merged or bundled envelope is not a newly computed nucleotide alignment.
Estimated lengths retain their flags; unobserved identity is never invented.

## Density is an explicit statistic {#density}
{{code:ribbondensity}}

{{figure:density}}

Each end is divided into four position bins. Ribbons darken when more alignments share the same sequence pair, direction and bins. The legend shows the computed count from `after_stat(bundle_n)`. Each input row still gives one ribbon. Use `stat_ribbon_bundle()` to aggregate neighbouring ribbons.

## Choose interval ribbons or point links {#points}

{{code:pointlinks}}

{{figure:pointlinks}}

These arrows join the midpoints of the first three alignments and point toward the subject. Use ribbons to show complete intervals.

## Check your result {#check-result}
Inspect each report and member rows before reading the plot. Aggregated bands are not new alignments; point links summarize midpoints of existing HSPs.

<div class="next-step">Continue: <a href="views.html">Focus, segments and homology →</a></div>
