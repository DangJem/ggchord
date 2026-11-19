Fusion combines plots that already have their own coordinate systems. Here a circle shows the full MT108731.1 genome, and two vertical backbones show it beside MT118296.1. An observed HSP will connect one interval on the circle to one on the second backbone.

<a id="draw"></a>

## Place two independently built plots {#draw}

Build the circular overview and linear detail first. Pass the overview to `ggchord_fusion()`, then place the second plot on the right with `fusion_panel()`. The central panel has the ID `center`. Give the new panel the ID `detail` so another layer can address it.

{{example:fusion-panels}}

<a id="adapt"></a>

## Route one HSP to the named panels {#adapt}

Use the first query-to-subject HSP again. In `fusion_ribbon()`, `qpanel` places MT108731.1 positions 26,298–51,139 bp on the circle. The `spanel` argument places MT118296.1 positions 7,121–31,959 bp on the right. Although MT108731.1 also appears in the right panel, the query end meets the circle.

{{example:first-fusion}}

<span id="shared-guides"></span>

Each panel must contain the named reference and the interval it receives. To collect matching outside legends across panels, use `guides = "collect"` in `ggchord_fusion()` and give comparable measures the same scale limits and breaks. Next, follow [the six-panel case](case-fusion.md) from a genome overview to a short alignment, or [export the figure at its final width](output.md).
