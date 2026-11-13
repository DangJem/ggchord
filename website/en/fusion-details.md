{{code:courseblockbr}}

{{prepare:fusiondata,fusionpanels,fusionlinear,fusionview}}

## A source interval and its local detail {#magnify}

`magnify_ggchord()` composes a Circular/Chord/Collinear source and an explicit Collinear or registered character detail plot. A source interval is biological data, not a rectangle in image pixels. The same region must be fully visible in both panels. Both rebuild at the final device size; no screenshot is stretched.

This example selects pETDuet-1 source positions 101–1900 (1800 bp). Its origin-crossing annotation at 5404→2 has no overlap with this interval; exclude it explicitly before linear focusing rather than taking its min/max as a long feature. The detail retains 10 annotations and labels the explicitly selected CDS/replication-origin subset. Reversing the detail's display orientation changes neither source coordinates nor biological strands.

{{code:courseblockbs}}

![A local source detail](../assets/coursefocusedview.png)

Use ordinary `ggsave()` at a chosen physical size, for example 12×8 inches; `view_ggchord()` uses the same renderer. `center_size` and `detail_size` are container fractions, while `highlight_width` is physical mm. Changing the detail side or orientation can change leader crossings; global collision avoidance is not claimed. `connector = "none"` keeps only highlights.

The dashed lines are decorative focus leaders. They are exported separately from biological `connections` and `relationship_data`; adding actual `fusion_ribbon()` layers still requires real relationship data. `focus_region_data` preserves the original selection and arbitrary metadata, while `focus_regions` reports resolved instances, source origins and the sampled physical path-length ratio. The ratio can be below one with a small detail allocation; it is not a change to biological length.

The detail plot owns its data, scales, layers and statistics. For count/bin/coverage tracks, prepare and state the intended local or global statistical population yourself; the wrapper records `explicit_detail_plot` and never silently chooses a scope. It does not copy arbitrary source layers, infer alignments or accept unregistered Cartesian character panels. Use one ordered region, explicitly select `instance`/`detail_instance` when an accession is repeated, and declare `wrap = TRUE` for origin-crossing selections as described below. Explicit nested/free-position allocations are described below; measured obstacle routing is described below.

## Paired source windows with real HSPs {#paired-detail}

In this example, each row of `region` selects one original source interval; all rows share one explicit detail panel. The example selects the two endpoints of a real HSP from the four-sequence dataset. Gene direction and the local HSP come from the focused input tables. `region_row` links each derived region, front and decorative leader to its input row; arbitrary original metadata stays in `focus_region_data`.

{{code:courseblockbt}}

![Both HSP ends in one detail](../assets/coursepaireddetail.png)

This does not create multiple insets or infer base-level alignments. Each interval must be fully visible in both panels. Instance selectors apply per row.

## Independent detail panels {#multiple-details}

Split the two endpoints of one real HSP into independent local views. A named
list supplies the detail plots, `detail_panel` routes each source row, and a
named `side` vector assigns distinct sides. Each panel owns its data and scales.

{{code:courseblockbu}}

![Two independent source detail panels](../assets/coursemultiple.png)

This example uses four named details, one per cardinal side; diagonal sides and explicit boxes are available below. Every panel needs a region; missing
or unknown targets fail. Source and detail instance selectors remain independent.
`detail_size` accepts a named list of width/height fractions; NULL uses Fusion
defaults. Overlaps fail rather than displacing explicit allocations. Standard
saving rebuilds panels and leaders at the requested size. At small sizes, inspect
labels and enlarge the canvas as needed. Leaders may cross the main plot or one
another; global avoidance is separate from the explicit nested/box allocation described below.

## Character details with exact source anchors {#character-detail}

A Cartesian character plot can be an independent detail panel after `register_ggchord_alignment()`. Register the actual displayed character rows: x is the original alignment/coding column, y is the alignment row ID, and label is the actual symbol. The source map remains separate from the ordinary ggplot scales. These real 60-column excerpts retain the upstream gaps and mismatches; their macro lengths refer to the **ungapped aligned fragments**, not complete genomes.

{{code:courseblockbv}}

![Registered details from actual aligned characters](../assets/coursecharacterdetail.png)

The ordinary plot remains editable: local colours, layers and annotations remain independent. Identity/reverse x scales and discrete y rows are supported. Each Fusion build verifies the registered residue positions/labels after edits and remeasures the final panel in mm; half-cell boundaries attach to the exact endpoint cells. Gap cells have no source position and are never interpolated into bases. A source interval must contain every selected base in biological order; split introns into separate source rows; origin joins use the explicit wrapping declaration below. Use `detail_instance` to select an alignment row when several rows share a source. Unregistered ggplots, facets, nonlinear/flipped coordinates and moved registered cells fail.

`translate_ggchord_cds()` nucleotide rows can be registered in the same way. Keep the translation on codon centres as in `examples/m2-detail.R`; use `source_positions` for full codon membership across segments. AA residue rows use their own aa coordinate domain; connecting a genomic bp macro directly to aa positions still requires an explicit CDS projection. Decorative magnify leaders remain separate from biological Fusion links. Character endpoints have no backbone colour for a Fusion gradient: supply explicit cross-panel fill. Narrow panels retain `geom_residue()`'s tile/error policy; enlarge the physical export instead of relying on DPI. Global leader avoidance is separate; explicit nested/free-position allocation is supported below.
