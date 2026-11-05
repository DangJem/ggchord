Density colouring and bundling are different operations. Retain original HSPs before declaring bins, weights, groups and aggregation meaning.

## Observed input and minimal fields {#input}

Independently use bundled complete sequence lengths and observed HSPs; no new alignment is computed. In addition to IDs/endpoints, weight=length requires length and weight=pident requires both length and pident. Do not infer those statistics from spans.

{{code:comparedensitydata}}


## Four versus eight bins: one ribbon per observation {#density}

Normalize each endpoint midpoint by its complete sequence length and bin by pair, direction and optional group_by. Changing only bins from 4 to 8 retains input HSP rows while changing grouping granularity. Bundle_n counts grouped observations, not sequencing depth; after_stat reads the computed result.

{{code:comparedensityfour}}
{{figure:comparedensityfour}}
{{code:comparedensityeight}}
{{figure:comparedensityeight}}


## Bundling combines displayed members {#bundle}

`stat_ribbon_bundle` creates envelopes within bins/groups with min_bundle at least 2. Smaller groups remain individual ribbons; input counts are not bundle counts. Envelopes are not newly computed nucleotide alignments and do not invent identity. `bundle_ggchord_ribbons` returns retained tables and reports for inspection before plotting.

{{code:comparebundle}}
{{figure:comparebundle}}
{{code:comparebundlesummary}}
{{output:comparebundlesummary}}


## Computed variables, weights and preprocessing {#parameters}

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



## Your alignments and exercises {#adapt}

Set bins to 6, observe colour changes and confirm retention of original observations. Preserve actual condition columns for group_by rather than pooling samples implicitly. Link statistics support chord/collinear views; single circular mode rejects link layers. Continue with [link preparation](ribbons.html), [real cases](real-data.html) and [scale families](scale-families.html).
