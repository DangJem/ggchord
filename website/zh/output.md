<a id="size"></a>

## 按最终宽度保存 {#size}

示例以 160 × 110 mm 保存四序列图，在当前 R 会话的临时目录写入 `sequence-backbones.pdf`。用于论文时，可将 `output_file` 改为长期保存路径。显式传入图对象，可以准确指定要保存哪一幅图。

{{example:export-size}}

<span id="layout-results"></span>

<a id="preview"></a>

## 按相同尺寸检查文字和间距 {#preview}

要静态预览，可给 `view_ggchord()` 传相同的图、宽、高和单位，`dpi` 控制栅格分辨率。最终用 `ggsave()` 保存，并按目标设备选择文件扩展名。印刷用 PNG 可显式设置 `dpi`，仍须保留要求的物理尺寸。

查看画布边缘的标签、面板标题与轴的距离，以及图例键是否清晰。字体度量与可用空间共同决定文字能否放下。需要留下测量结果时，见[布局诊断](diagnostics.md)。
