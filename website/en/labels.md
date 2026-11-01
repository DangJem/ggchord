<div class="task-intro"><p>Arrange text while retaining annotation geometry: choose export dimensions and a label subset before choosing automatic or fixed placement.</p></div>

## Prepare this page {#prepare}
In a fresh R session, run the code below, then follow this page in order. Colours and layout are optional display choices; inputs select only the fields needed.

{{prepare:setup,docstyle}}

## Select a layout {#layout}
Gene labels offer `"auto"`, `"radial"` and `"arc"` layouts. General feature
labels additionally support the internal-first `"feature"` default and external
`"callout"` layout. Choose the documented values for the layer you use; they are
not interchangeable presets.

After this page’s preparation code:

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

## Match anchors across layers {#anchors}
In collinear mode, an omitted feature-label Position inherits an unambiguous
matching feature layer. Unique `feature_id` permits reordering or subsets, but
accession, interval and strand must agree. Multiple matching layers need an
explicit Position or unambiguous IDs. Preview after selecting the final device
size; resizing the Plots pane changes the available space.

## Compare the three gene layouts {#comparison}
{{code:labellayouts}}
{{figure:labelauto}}
{{figure:labelradial}}
{{figure:labelarc}}

Keep at most the first and last informative labels per sequence, retaining all gene arrows, to compare text orientation.
Print `label_plots$auto`, `$radial` or `$arc` to inspect a candidate. The inputs and
requested font size are identical; only the layout choice changes. Auto uses
shared columns where applicable, radial follows local contours and arc follows
the sequence tangent. `p_manual_labels` shows the separate fixed-label route;
its explicit rotation does not belong in an automatic repel layer.
General feature labels have additional internal-first/callout choices; see the
[complete circular workflow](plasmids.html) for dense mixed annotation and
[stacking](annotations.html#stack) for inherited or explicit Position anchors.

## Check your result {#check-result}
The three layout comparisons share the same input. Assess text direction, leaders and hiding reasons together; missing text does not mean a feature was lost.

<div class="next-step">Continue: <a href="output.html">Export, preview and diagnostics →</a></div>
