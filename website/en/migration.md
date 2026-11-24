<span id="api"></span>

Use v0.29.0 for scripts written against the old API. From v0.50.0 onward, ggchord uses a redesigned API that is incompatible with the old interface.

For scripts written against earlier development snapshots of v0.50.0, put `layout_label_repel()` in `label_layout` on your usual label layer. Use `layout_label_fixed()` for fixed text and `geom_seq_label_center()` for circle-centre text. Move `leader` and its `segment_*` styles to the label layer. Set sequence-axis visibility, gap and text orientation through `scale_seq_continuous()` and `guide_seq_axis()`.

`coord_circular()` now measures `start_angle` counterclockwise from twelve o’clock. Subtract 90 degrees from an earlier v0.50.0 development angle to preserve its orientation. Use `start_angle = -90` to start at three o’clock.

Use `"arrow_arc"` wherever an earlier development snapshot used `"primer_arc"`, including fixed annotations, `shape_feature()` and manual shape scales. The `geom_primer()` shortcut uses the renamed preset by default.
