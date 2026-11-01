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

## 组合柱形、折线与类别 tile {#multiple}
{{code:multitrack}}
{{figure:multitrack}}

柱形和折线有意共享 `signal` 及一个数值域；折线是合成信号的一半，不是另一份实测
样本。两层声明相同宽度。类别轨道 `sites` 有独立 ID 和填色标度，零基线参与数值域
训练。显式 limits 应表达目标测量范围，调整数值 scale 不移动轨道空间位置。
`focused_bins` 为选定 locus 保留区间裁剪和源映射。检查轨道诊断中的边界裁剪，
不能把遗漏行当作零。轨道颜色与填色角色见[样式](styling.html#roles)。

折线位置使用向下取整的整数 bin 中点。第二条分类带明确位于外侧，排在数值带之后；
小半径下过宽的内侧轨道可能折叠，此时会报错，不会穿过圆心继续绘制。
