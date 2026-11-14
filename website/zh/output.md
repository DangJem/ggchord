## 本页准备 {#prepare}
用四条随包序列画出共线性骨架，后面用同一张图演示预览和保存。

{{code:frontlinearbase}}
{{figure:frontlinearbase}}

## 按目标尺寸保存 {#size}
下面以本页创建的 `p_linear_base` 为例。有自己的图时，把它换成你的绘图对象。

{{code:frontpreview}}

圆图省略宽高可自动估计尺寸。弦图和共线性图可指定宽度并设 `height = NULL`，按图形估计高度。网页缩略图还会按正文栏缩放；仅增大 dpi 不会让网页中的字变大。

确认字号、标签和图例后，使用相同尺寸保存。

{{code:frontsave}}

`view_ggchord()` 返回临时图片路径，不保存最终导出文件。
网站、README 和 vignette 共用由 `view_ggchord()` 生成的预览 PNG，尺寸见[图集](gallery.html)，
不截取临时调整后的交互面板。最终矢量输出应显式选择目标 PDF/SVG 设备；字体可用性
和嵌入仍取决于设备。

## 检查要点 {#check-result}
预览与保存使用相同宽高。检查文字、图例与警告；一次测量成功不代表跨设备字体和发布质量已认证。

## 选择文件格式 {#formats}

PNG 适合网页与幻灯片；PDF/SVG 用于需要矢量输出的工作流。文件名扩展名选择默认设备，明确指定设备时以该设备为准。宽高继续复用前面的目标尺寸；字体可用性与嵌入须在最终文件中核对。

```r
ggsave("comparison.pdf", p_linear_base,
  width = output_width, height = output_height, units = "in")
```
