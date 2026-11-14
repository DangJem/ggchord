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

Choose point versus interval geometry in the [links section](multi-sequence-links.html#interval-points).

## Result checks {#check-result}
Inspect each report and member rows before reading the plot. Aggregated bands are not new alignments; point links summarize midpoints of existing HSPs.

<div class="next-step">Continue: <a href="views.html">Focus, segments and homology →</a></div>

Continue to [link statistics](ribbons.html#statistics) for bundle_n, bundle_weight, density and bins/weight/group_by/min_bundle.

Density colouring and bundling are different operations. Retain original HSPs before declaring bins, weights, groups and aggregation meaning.

## Observed input and minimal fields {#statistics}

Independently use bundled complete sequence lengths and observed HSPs; no new alignment is computed. In addition to IDs/endpoints, weight=length requires length and weight=pident requires both length and pident. Do not infer those statistics from spans.

{{code:comparedensitydata}}


## Four versus eight bins: one ribbon per observation {#statistics-density}

Normalize each endpoint midpoint by its complete sequence length and bin by pair, direction and optional group_by. Changing only bins from 4 to 8 retains input HSP rows while changing grouping granularity. Bundle_n counts grouped observations, not sequencing depth; after_stat reads the computed result.

{{code:comparedensityfour}}
{{figure:comparedensityfour}}
{{code:comparedensityeight}}
{{figure:comparedensityeight}}


## Bundling combines displayed members {#statistics-bundle}

`stat_ribbon_bundle` creates envelopes within bins/groups with min_bundle at least 2. Smaller groups remain individual ribbons; input counts are not bundle counts. Envelopes are not newly computed nucleotide alignments and do not invent identity. `bundle_ggchord_ribbons` returns retained tables and reports for inspection before plotting.

{{code:comparebundle}}
{{figure:comparebundle}}
{{code:comparebundlesummary}}
{{output:comparebundlesummary}}


## Computed variables, weights and preprocessing {#statistics-parameters}

| Parameter/variable | Default or meaning |
| --- | --- |
| bins | 80; positive integer normalized midpoint bins per sequence |
| weight | length; length/count/pident determines weight source |
| group_by | NULL; columns defining independent conditions |
| min_bundle | 2; bundling requires at least two members |
| bundle_n | Original member count within a group |
| bundle_weight | Sum of row weights: count uses 1, length uses observed length, pident uses length × pident / 100 |
| density | Group weight divided by the maximum group weight for that ordered sequence pair, in 0–1; not absolute coverage depth |

Filtering selects observations by criteria; deduplication selects representatives; merge combines compatible neighbourhoods; bundle aggregates normalized bins. Order changes the input set, so retain reports and membership at each step. Envelopes are not original alignments for unqualified biological inference.



## Use your alignments {#statistics-adapt}

Set bins to 6, observe colour changes and confirm retention of original observations. Preserve actual condition columns for group_by rather than pooling samples implicitly. Link statistics support chord/collinear views; single circular mode rejects link layers. Continue with [link preparation](ribbons.html), [real cases](real-data.html) and [scale families](scale-families.html).
