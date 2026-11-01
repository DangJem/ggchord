## 声明位置与数值域 {#domain}
先运行[共享准备代码](index.html#quick-start)：

{{code:tracks}}
{{figure:tracks}}

信号是任意单位的合成正弦波，不是测量得到的覆盖度或 GC skew。
ggchord 绘制已提供的数值，不计算序列组成。替换成实际区间值后应声明真实单位。

## 有意共享轨道 ID {#sharing}
`position_track("signal", ...)` 选择位置与物理宽度，
`scale_track_value_continuous("signal", ...)` 定义数值域，不负责空间放置。
同 ID 图层联合训练一个数值域，包括柱形基线，且轨道声明必须一致。
填色/颜色 scale 另按该 ID 限定作用域。`geom_track_tile()` 表达类别，
`geom_track_line()` 连接数值观测，`geom_track_bar()` 表达相对于基线的区间值。

## 检查保留结果 {#diagnostics}
通过 `export_ggchord_layout()` 查看 `track_registry`、`track_diagnostics`。
原始数值范围、物理 mm 贡献与原生位移的单位各不相同。
使用 `focus_ggchord_track()` 按相同 loci 裁剪轨道，检查裁剪或遗漏行，不能把未显示的柱形当成零。
