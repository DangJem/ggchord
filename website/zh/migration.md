<span id="api"></span>

使用旧 API 的脚本请选用 v0.29.0。从 v0.50.0 起，ggchord 采用重构后的新 API，与旧接口完全不兼容。

对于 v0.50.0 早期开发快照编写的脚本，请在原标签层的 `label_layout` 中传入 `layout_label_repel()`。固定文字使用 `layout_label_fixed()`，圆心文字使用 `geom_seq_label_center()`。将 `leader` 及 `segment_*` 样式移到标签层。通过 `scale_seq_continuous()` 和 `guide_seq_axis()` 控制序列轴的可见性、间距及文字方向。

`coord_circular()` 现在从十二点方向逆时针量取 `start_angle`。将 v0.50.0 早期开发快照的角度减去 90°可保留原朝向。要从三点方向起画，设 `start_angle = -90`。

如果早期开发快照在固定注释、`shape_feature()` 或手动形状标度中使用 `"primer_arc"`，请改用 `"arrow_arc"`。`geom_primer()` 快捷层默认使用改名后的预设。
