Fusion combines plots that already have their own coordinate systems. Here a circle shows the full MT108731.1 genome, and two vertical backbones show it beside MT118296.1. An observed HSP will connect one interval on the circle to one on the second backbone.

<a id="draw"></a>

## Place two independently built plots {#draw}

Build the circular overview and linear detail first. Pass the named list `list(overview = overview, detail = detail)` to `ggchord_fusion()`. The default grid arranges the plots in that order. The names are the panel paths used by layers that connect them.

{{example:fusion-panels}}

<a id="adapt"></a>

## Route one HSP to the named panels {#adapt}

Use the first query-to-subject HSP again. Add `qpanel_path` and `spanel_path` to the relationship data for `fusion_ribbon()`. The query goes to MT108731.1 at 26,298–51,139 bp in `overview`, and the subject goes to MT118296.1 at 7,121–31,959 bp in `detail`. The accessions identify the source sequences and the panel paths select the intended views.

{{example:first-fusion}}

<span id="shared-guides"></span>

Each panel must contain the named reference and the interval it receives. To collect matching outside legends across panels, use `guides = "collect"` in `ggchord_fusion()` and give comparable measures the same scale limits and breaks. Next, follow [the six-panel case](case-fusion.md) from a genome overview to a short alignment, or [export the figure at its final width](output.md).
