Fusion places independently constructed plots on one canvas. A cross-panel ribbon connects its two original sequence intervals to named panels, letting an overview and a detailed view use different coordinates.

<a id="draw"></a>

## Connect a circle to a linear view {#draw}

The centre plot contains MT108731.1. The right panel shows both sequences. Their shortest observed alignment connects 20,376–20,677 bp on MT108731.1 to 403–704 bp on MT118296.1. `qpanel` places the query end on the circle and `spanel` places the subject end in the panel called `detail`.

{{example:first-fusion}}

<span id="shared-guides"></span>

<a id="adapt"></a>

## Choose an unambiguous destination {#adapt}

Check that both endpoints fall inside the visible interval of the intended panel and use compatible units. The side chosen in `fusion_panel()` places the child plot. Its coordinate system controls the direction of the child backbones. For an overview connected to a source window, the magnification helper can create and register the window.

To place outside legends together, pass `guides = "collect"` to `ggchord_fusion()`. Give panels matching limits, breaks and styles when they represent the same quantity. The composition collects their legends while retaining each panel’s scales.

[Six-panel comparative map](case-fusion.md) · [Cross-panel ribbons](../../man/fusion_ribbon.Rd) · [Magnified source windows](../../man/magnify_ggchord.Rd)
