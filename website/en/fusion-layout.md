{{code:courseblockbm}}

{{prepare:fusiondata,fusionpanels,fusionlinear,fusionview}}

## Styles and physical layout {#styles}

Cross-panel layers accept `aes(fill=..., colour=..., alpha=..., linewidth=...,
linetype=...)` with ordinary ggplot2 scales. Panel role-specific scales remain
local. Explicit fixed styles override their mappings. Default Identity is
0–100 with 0/50/80/90/95/100 proportional breaks. Endpoint gradients and
subject tips are opt-in, using the same renderer as single-view ribbons.
A tip is a decorative direction marker; full source fronts remain exported.

{{code:fusionstyle}}

Default `layout="fixed"` allocates panel rectangles and rebuilds their labels
and guides at the output size. `layout="adaptive"` additionally fits the
central rectangle to the physical aspect ratio and jointly budgets native outer Guides. Explicit font sizes remain
fixed. This is not a general font-scaling or global collision solver. Each
panel retains its own guides; different scale meanings are not automatically
merged. Use `show.legend=FALSE` or local `theme()` deliberately where appropriate.
Composition `guides()` configures cross-panel Guides; `theme()` also styles the native outer plot furniture described below. Required
guide space is measured; a canvas too small for it fails with a diagnostic.

## Eight directions and independent physical placement {#advanced-placement}

`fusion_panel()` supports the four cardinal sides plus `topright`, `topleft`, `bottomright` and `bottomleft`. Default corner layouts reserve space for all eight; side names remain unique. For other arrangements use `angle`/`distance` or an explicit `box`. There is no fixed panel-count limit for explicit boxes/angles: every allocation must still fit. Placement changes the container, not the child's Coord or text orientation. Angles describe physical directions on rectangular exports; distance is a fraction of half the shorter allocation dimension, measured from the centre allocation.

{{code:courseblockbn}}

![Explicit layout in millimeters](../assets/coursephysicalbox.png)

A numeric box uses xmin, xmax, ymin, ymax fractions; an absolute-unit box keeps its physical dimensions. Coordinates refer to the composition content area after native outer plot decorations and cross-panel Guides. Do not combine box with side/size/angle/distance, or center_box with center_size. Physical allocations are checked at every target size; they are never silently clamped. Explicit `allow_overlap = TRUE` on `ggchord_fusion()` or `magnify_ggchord()` permits deliberate overlays. Later panels cover earlier panels, using vector grid clipping rather than PNG crops; translucent connections are not repainted twice in the same area. The export's `overlaps` table records covering order. This permits intentional masks and does not promise automatic obstacle avoidance.

Sankey source IDs occur in both mirrored columns. Choose an exact `q_instance`/`s_instance` from `get_chord_layout(plot)$display_instances`; omission fails when ambiguous. Radial and Chord use their own registered backbone frames. All q/s intervals, units and source records remain unchanged. If a declared `source_hsp_id` appears repeatedly, every participating Fusion layer must opt into `allow_repeated = TRUE`, including layers inside nested compositions. `source_observations` distinguishes explicit evidence IDs from display counts. This does not infer identity by comparing equal table values. Dotplot has independent axes and is not yet a registered Fusion endpoint.

## Outer Guides and the panel budget {#outer-guides}

Keep each panel's local scales and Guides independent. Composition scales and
`guides()` describe the cross-panel ribbons. Native Guides may occupy different
outside positions, rather than all requiring one shared side:

{{code:courseblockbo}}

![Multi-side guide allocation](../assets/courseguidefusion.png)

`layout = "fixed"` keeps existing placements. `"adaptive"` reuses the native
Guide candidate/measurement engine, jointly considers side budgets and panel
allocations, and then rebuilds every child in the remaining physical space.
It may move/reflow unlocked Guides; explicit position, direction, rows/columns,
fonts, key dimensions and absolute panel boxes remain constraints. For example,
a default Identity ruler can become horizontal on a wide, short canvas; its
length is remeasured for all six original labels rather than dropping 95 or 100.

`metadata$guide_layout` reports choices and locks. `guide_bounds_mm` and
`guide_sizes_mm` describe the final outer boxes, after packaging and before
panel measurement. Nested composition Guide reports appear in `diagnostics`.
Impossible constraints fail explicitly. `theme(legend.position = "none")`
hides composition Guides using native ggplot2 semantics. Inside composition
Guides remain unsupported; each independent panel can still own local Guides.
This is a bounded physical-space search, not a claim of global optimality.

## Automatic perimeter panels {#auto-perimeter}

Set `fusion_panel(..., side = "auto")` to place more than eight independent
panels around the center. The placement uses actual mm after the outer Guide
budget, clockwise in input order on rectangular canvases. Explicit sizes and
all fixed named/angle/box allocations stay fixed. Unlocked default containers
may adapt down to a 15 mm floor; this does not shrink fonts or guarantee that
arbitrary text will fit. Use each child's own label/physical-size policy.

{{code:courseblockbp}}

![Twelve deliberate views of one HSP](../assets/courseautofusion.png)

This real example intentionally displays one observed HSP twelve times. Its
explicit `source_hsp_id` and `allow_repeated = TRUE` distinguish 12 connections
from one observation. Each view still rebuilds its exact interval fronts.
`metadata$panel_layout` records the bounded search, preserved constraints and
container adaptation; panel rows expose `auto_angle` and `auto_size_factor`.
Multiple magnify details can use repeated named `"auto"` sides, with explicit
`detail_panel` selections. If no packing is feasible, the build fails with a clear diagnostic;
no panel is discarded. The search is deterministic and conservative, not a
global optimum or a promise that any count fits every canvas.

## Native outer plot furniture {#plot-furniture}

Use ordinary `labs()` and plot Theme elements for the composition's title,
subtitle, caption, tag, background and margins. These settings previously
entered the style carrier without being drawn; they now reserve measured
space and appear in the final static drawing. Each child retains its own
axes, theme and source coordinates. Fusion has no Cartesian data axes.

{{code:courseblockbq}}

![Outer title and caption](../assets/coursestyledguidefusion.png)

`metadata$content_bounds_mm`/`content_canvas_mm` describe the actual inner
composition frame; `plot_decoration` records effective labels and native
physical envelopes. Boxes refer to this content frame after outer decoration
and Guide allocation. Existing fixed boxes are rechecked at target size;
enlarge the canvas or revise explicit margins if they no longer fit. Clear
headings with native `labs(title = NULL, subtitle = NULL, caption = NULL, tag = NULL)`.
No decoration keeps the previous rendering path.

Nested compositions retain independent headings and mm registration. Native
foreground/background separation and inherited clipping keep translucent
connections from accumulating repeated paints in explicit overlap mode.
The source tables, selected intervals and palette mappings remain unchanged.

Ordinary child margins do not act as masks for composition connections. In explicit overlap mode, links traverse the full container below its foreground; later containers retain the declared covering order. Use measured `route = "avoid"` for annotation clearance instead of relying on accidental margin gaps.

## Integrated composition: global views and a local detail {#flagship}

Four independently valid views retain common sources. `box` controls physical allocation; `focus_ggchord_data()` explicitly prepares local windows. Cross-panel ribbons reuse source HSPs, and every participating layer declares repeated display.

{{code:flagshipdata,fusionoverview,flagshipfusion}}

{{manualfigure:flagfusion}}
