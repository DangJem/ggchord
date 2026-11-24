Fusion combines plots that already have their own coordinate systems. Here a circle shows the full MT108731.1 genome, and two vertical backbones show it beside MT118296.1. An observed HSP will connect one interval on the circle to one on the second backbone.

<a id="draw"></a>

## Place two independently built plots {#draw}

Build the circular overview and linear detail first. Pass the named list `list(overview = overview, detail = detail)` to `ggchord_fusion()`. Add `layout_fusion_grid(ncol = 2, gap = 3)` to arrange the views in one row with a 3 mm gap. The names are the panel paths used by layers that connect them.

{{example:fusion-panels}}

<a id="adapt"></a>

## Route one HSP to the named panels {#adapt}

Use the first query-to-subject HSP again. Map `qpanel_path` and `spanel_path` to the two view names in `geom_link_ribbon()`. Use `geom_link_line()` for observed point-to-point relationships. The query goes to MT108731.1 at 26,298–51,139 bp in `overview`, and the subject goes to MT118296.1 at 7,121–31,959 bp in `detail`. The accessions identify the source sequences and the panel paths select the intended views.

{{example:first-fusion}}

<span id="shared-guides"></span>

Each view must contain the specified source positions. Use `fusion_magnify_auto()` to derive a window from a sequence overview or `fusion_magnify_manual()` to mark the same source region in two prepared views. With real characters, Auto can show a residue detail when the panel fits the measured glyphs. Set `level = "feature"` to keep the annotation window or `level = "residue"` to request a character detail. For the latter, retain the reference strings with `read_fasta_lengths(keep_sequence = TRUE)` or start from a registered character plot. Set `guides = "collect"` in `ggchord_fusion()` to arrange outside legends together. Use `scale_for()` with a table of `panel_path` and `layer_id` values to train a visual scale across selected views. Next, follow [the six-panel case](case-fusion.md) from a genome overview to a short alignment, or [export the figure at its final width](output.md).
