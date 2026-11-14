## Locate the problem by symptom {#symptoms}

Check the input and current plot before inspecting diagnostics. These actions preserve original observations and coordinates.

| Symptom | Minimal check | Adjustment and verification |
| --- | --- | --- |
| Missing sequence or feature | Exact accver matching, valid complete lengths and units | Correct evidence-backed ID mappings and rebuild a minimal plot; see [table checks](data.html#validation) |
| Missing link | Endpoint IDs, intervals, explicit filters and the Coord's adjacency default | Use `adjacent_only = FALSE` when all-pair selection is intended; inspect [selection rules](multi-sequence-views.html#selection) and verify the source row in exports |
| Crowded or hidden labels | Run `diagnose_ggchord()` at final dimensions; inspect visibility reasons and the plot | Adjust [label layout](labels.html), wrapping or canvas; recheck visibility/conflicts without treating missing counts as zero |
| Exported positions differ or guides clip | Compare preview/save width, height, units, device, fonts and warnings | Reuse [target dimensions](output.html#size) and inspect the saved file; dpi does not replace physical dimensions |

Requires ggchord and optional dplyr for later tables. Start with a minimal backbone at a normal size.

{{code:courseblockbj}}

{{prepare:frontlinearbase}}

```r
report <- diagnose_ggchord(p_linear_base, width = 8, height = 5, device = "png")
report$checks
```

Inspect failures/warnings and use the matching action above. Later sections explain export fields and measurement boundaries as needed.


## Consume geometry by name {#geometry}
{{code:frontgeometry}}

{{output:frontgeometry}}

`export_ggchord_layout()` is the recommended layer-aware table interface.
`get_chord_layout()` retains the complete documented layout; neither is deprecated.
Select named fields, not fixed column positions. `source_row` indexes resolved
layer input; IDs are build-local resources, not permanent biological identifiers.
See the complete [layout contract](reference/ggchord_layout_contract.html).

## Interpret evidence {#diagnostics}
`diagnose_ggchord()` separates measured panel size from geometry convergence.
A measured panel does not certify collision-free text, font embedding or all
applications. Inspect warnings, text/guide measurements and placement reasons.
Drawing warns when bounded physical refinement does not converge. Geometry
measurements then remain estimates; diagnostics capture the warning and retain
the `physical_layout` check. Try a larger canvas or revise label/guide settings,
then measure again.
`optimize_ggchord_layout()` applies documented opt-in adjustments; retain a copy
and inspect the result. Previews and diagnostics use independent render copies.
Direct mutation of shared ggproto layer environments retains ordinary ggplot2
reference semantics; it is not an automatically isolated edit.

## Complete layout, tables and opt-in ordering {#consumption}

The preceding table identifies the backbone layer. Set `original_data = TRUE` in `export_ggchord_layout()` to include input rows. See the [layout contract](reference/ggchord_layout_contract.html) for the full layout and [order and orientation](reference/optimize_ggchord_layout.html) for reordering.

## Measure label contacts {#label-quality}

`diagnose_ggchord(p, width = 10, height = 8)$label_quality` reports measured panel text objects, annotation leaders, box contacts and physical leader lengths. `summary` distinguishes source labels retained/hidden by layer policy from text objects and glyph boxes. `conflicts`, `text` and `leaders` locate each reported contact. Its `measurement` entry keeps panel measurement, geometry convergence and physical refinement passes separate. Circular refinement can test one bounded scale prediction after two measurements; `refinement_acceleration` records the outcome. A prediction is accepted only after rebuilding and passing the existing measured convergence gate. Rejection preserves the ordinary sequence and does not become a pass. Ordinary refinement uses a shared bounded budget across views. Glyphs within one arc string share a text-object identity; touching polyline endpoints are not counted as crossings. Inspect actual PNG output as well as the tables.

Registered text uses the active target device for physical measurement, including ordinary labels and arc-glyph advances; a layout without an active device retains its measurement fallback. Font fallback and embedding still need separate review. Restriction labels use the same bold-name/plain-coordinate fragments and physical gap for drawing, perimeter packing and measured routing. Their effective theme or explicit type settings enter layout before placement. These are oriented text envelopes, not exact ink outlines. A leader crossing a box can need visual review; whitespace and label backgrounds are not certified glyph shapes. `feature_contacts` additionally checks registered gene/feature polygon envelopes, including transparent ones. It separates proven own-source label placement (`own_source`), contacts with other constructor records (`other_source`), and associations without sufficient evidence (`unresolved`). Source identity is never guessed from equal contents. `feature_polygons` provides the measured vertices; `feature_parts_unchecked` lists skipped non-polygon strokes/markers. Stroke width, visibility/occlusion, outer decorations and unregistered grobs remain unchecked; a zero polygon-contact count does not certify them. Diagnostics do not improve a layout by hiding labels or changing source data. A measured panel or zero contacts in the checked categories does not establish that all layout constraints passed.

SnapGene serves as a visual reference for clear annotation hierarchy and orderly leaders; matched exports from that application are not required. Use controlled ggchord comparisons to assess changes in content retention, physical size, crowding and legibility.

Automatic sequence-axis lines and major/minor ticks draw before the first gene, feature or primer layer. Opaque feature fills cover intersecting marks; a transparent feature remains transparent. Axis text keeps its measured clearance and draws in front. Your layers keep their relative addition order, including explicitly added `geom_seq_scale()` bars; tick coordinates, lengths and source features do not change.

With the usual overlap policy enabled, automatically positioned circular parallel labels, including seam strings such as `1 / 2856`, can gain bounded inward or outward radial clearance if their text envelopes touch a feature polygon or a fixed annotation-leader root. Other text and backbone geometry remain protected; internal and exterior annotation leaders then route around the resolved axis text. Mutable leader shoulders do not feed back into axis placement. Automatic endpoint strings may use either readable side of the same tick, with a tangential gap bounded at 4 mm, to avoid fixed feature exit rays. The `axis_clearance` table records the status, signed radial shift, `label_side` and `tangent_gap_mm`. Tick positions and label contents stay fixed; an explicitly supplied `axis.text.offset`, including a value equal to the preset, prevents this adjustment, as do explicit scale breaks. This is a local clearance rule, not a general axis-layout guarantee. Circular exterior feature and restriction-site leaders can then jointly try straight or single-elbow exits around the registered text. Source roots and label entries stay fixed; only a complete group with no text contacts or new leader crossings is accepted. Failed groups retain their existing routes and overlap appearance. Explicit straight and elbow restriction styles retain their selected routes.

Circular exterior callouts now protect the measured origin-axis and inner feature text on their way out of the feature band. Only paths with an oriented-box contact are rerouted; the true root, final label attachment and label contents stay fixed. The bounded search verifies its candidate against all registered static text and reports `leader_text_status` as clear, detoured or unresolved. An unresolved route retains its evidence for review. Unused zero-length placeholders are removed rather than treated as real leaders or marker observations.

On cramped circular canvases, already-crossed exterior leaders can acquire bounded detours around measured text and other leaders. Their real roots, text-edge attachments and label positions are retained. Unresolved routes remain visible in diagnostics; zero measured crossings is not a promise about every input or exact ink.

The `scale` and `adaptive` policies also adjust omitted fixed linewidths in native sequence, point-link, ribbon and alignment layers. Explicit linewidth values, including values equal to the default, remain fixed; mapped linewidths keep their scales. `fixed` and plots without a policy retain their styles. This does not multiply every geometry constant or native theme-derived aesthetic expression.

{{prepare:comparescaledata}}

## Target-size diagnostics and explicit ordering {#target-size}

`diagnose_ggchord` inspects one PNG size and returns measurement, text and guide hints without mutating the plot. `export_ggchord_layout` returns layer tables and original_data=TRUE retains source rows. `get_chord_layout` exposes the complete public layout. `optimize_ggchord_layout` returns proposed order/directions and before/after scores; it neither edits inputs nor applies changes automatically. Readability optimization is not biological correctness.

{{code:comparediagnostics}}
{{console:comparediagnostics}}
