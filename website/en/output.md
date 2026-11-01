## Save the intended size {#size}
After creating `p_collinear` in [three modes](modes.html):

{{code:save}}

RStudio's Plots pane uses the current panel size and can rearrange labels and
guides. Preview at the same explicit dimensions before final `ggsave()`:

{{code:preview}}

`view_ggchord()` returns a temporary image path; it does not keep a final export.
The website, README and vignettes share PNG assets exported through standard
`ggsave()` at sizes listed in the [gallery](gallery.html). They do not capture an
interactively resized panel. Use your intended PDF/SVG device explicitly for
final vector output; font availability and embedding remain device-dependent.

## Consume geometry by name {#geometry}
{{code:provenance}}

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
{{code:layoutresults}}

Create `p_collinear` in [three modes](modes.html#collinear) first.
`full_layout` retains the complete documented getter result. `layer_tables`
selects sequence/gene tables and resolved original data; it is not a raw-file
archive. Layer IDs and source-row roles distinguish independent same-kind layers.
Empty components can lack a fixed column prototype. Select fields by name.

`layout_choice` optimizes declared order/orientation using count weighting because
the minimal plotting table deliberately omits alignment-length statistics. It returns a choice
for explicit use in sequence settings; it does not mutate the plot or establish
biological orientation. Save the choice and compare the result before adopting it.
The export [asset manifest](../assets/plots.csv) records exact documentation sizes.
