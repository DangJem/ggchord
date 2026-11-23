<span id="api"></span>

Use v0.29.0 for scripts written against the old API. From v0.50.0 onward, ggchord uses a redesigned API that is incompatible with the old interface.

For scripts written against earlier development snapshots of v0.50.0, put `layout_label_repel()` in `label_layout` on your usual label layer. Use `layout_label_fixed()` for fixed text and `geom_seq_label_center()` for circle-centre text. Move `leader` and its `segment_*` styles to the label layer. Set sequence-axis visibility, gap and text orientation through `scale_seq_continuous()` and `guide_seq_axis()`.
