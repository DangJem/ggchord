<div class="task-intro"><p>Inspect figures at their actual export dimensions, save with standard ggsave(), and understand what layout tables and diagnostics establish.</p></div>

## Prepare this page {#prepare}
In a fresh R session, run the code below, then follow this page in order. Colours and layout are optional display choices; inputs select only the fields needed.

{{code:frontlinearbase}}
{{figure:frontlinearbase}}

## Save the intended size {#size}
The example uses `p_linear_base` from this page. If you already have a plot, use your plot object in its place.

{{code:frontpreview}}

For circular plots, omit both dimensions to estimate the preview size. For chord and collinear plots, supply a width and use `height = NULL` to estimate height. A webpage also scales a thumbnail to its content column, so increasing dpi alone does not enlarge its text.

After checking text, labels and legends, save at the same dimensions.

{{code:frontsave}}

`view_ggchord()` returns a temporary image path; it does not keep a final export.
The website, README and vignettes share PNG previews generated with
`view_ggchord()` at sizes listed in the [gallery](gallery.html). They do not capture an
interactively resized panel. Use your intended PDF/SVG device explicitly for
final vector output; font availability and embedding remain device-dependent.

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
`optimize_ggchord_layout()` applies documented opt-in adjustments; retain a copy
and inspect the result. Previews and diagnostics use independent render copies.
Direct mutation of shared ggproto layer environments retains ordinary ggplot2
reference semantics; it is not an automatically isolated edit.

## Complete layout, tables and opt-in ordering {#consumption}

The preceding table identifies the backbone layer. Set `original_data = TRUE` in `export_ggchord_layout()` to include input rows. See the [layout contract](reference/ggchord_layout_contract.html) for the full layout and [order and orientation](reference/optimize_ggchord_layout.html) for reordering.

## Check your result {#check-result}
Use matching dimensions for preview and saving. Inspect text, guides and warnings; a successful measurement does not certify fonts or release quality across devices.

<div class="next-step">Continue: <a href="faq.html">Frequently asked questions →</a></div>
