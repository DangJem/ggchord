<div class="task-intro"><p>Arrange text while retaining annotation geometry: choose export dimensions and a label subset before choosing automatic or fixed placement.</p></div>

<span id="prepare"></span>

## Choose a chord label layout {#layout}
Gene labels offer `"auto"`, `"radial"` and `"arc"` layouts. General feature
labels additionally support the internal-first `"feature"` default and external
`"callout"` layout. Choose the documented values for the layer you use; they are
not interchangeable presets.


{{code:labels}}
{{figure:labels}}

## Fixed and automatic placement {#fixed}
Use `geom_gene_label()` or `geom_feature_label()` for manual `label_rotation`,
radial/circumferential offsets and orientation controls. Automatic layers reject
these manual controls: choose `label_layout` or remove the manual argument.
Use `size`, `colour` and `family` for supported fixed text styles.

`label_fit`, `label_wrap`, `label_max_lines` and `label_min_size` control fitting
where supported. Shrinking is explicit and bounded. Unresolved labels may move
or hide; the feature geometry remains. Inspect `label_fit_diagnostics` and
`label_diagnostics` at the final export size.

{{code:manualgenelabels}}
{{figure:manualgenelabels}}

This fixes text rotation at 30 degrees. Fixed placement does not spread neighbouring labels; unresolved text may hide. One label is visible here while all four CDS arrows remain. Compare the automatic layouts below when more names need to be displayed.

## Match labels to arrows {#anchors}
Use the same CDS records for arrows and labels, preserving identifiers, intervals and biological strands. Filter labels from that table without changing source coordinates. Fixed layouts allow radial offsets and text rotation.

## Compare the three gene layouts {#comparison}
{{code:labellayouts}}
{{figure:labelauto}}
{{figure:labelradial}}
{{figure:labelarc}}

Compare the same four named CDS in each layout.
Pass `label_plots$radial` or `label_plots$arc` to `view_ggchord()` to inspect another layout at the same size. The inputs and
requested font size are identical; only the layout choice changes. Auto uses
shared columns where applicable, radial follows local contours and arc follows
the sequence tangent. `p_manual_labels` shows the separate fixed-label route;
its explicit rotation does not belong in an automatic repel layer.
General feature labels have additional internal-first/callout choices; see the
[complete circular workflow](plasmids.html) for dense mixed annotation and
[stacking](annotations.html#stack) for inherited or explicit Position anchors.

## Result checks {#check-result}
The three layout comparisons share the same input. Assess text direction, leaders and hiding reasons together; missing text does not mean a feature was lost.

<div class="next-step">Continue: <a href="parameters.html">Flexible parameter formats →</a></div>
