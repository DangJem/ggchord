先运行一次[共享准备代码](index.html#quick-start)。下列对象都可按普通 ggplot 组合，
打印或使用 `ggsave()` 保存。

## 弦图：比较多条序列 {#chord}
{{code:chord}}
{{figure:chord}}
`coord_chord()` 控制全局旋转和拟合。逐序列半径、曲率、间隙和方向属于 `geom_seq()`。
连接控制点为 NULL 时，弦图保留历史的原点路径，并不表示各模式统一使用直线。
ribbon 默认避障为 `uniform`，line 默认避障为 `none`。

## 圆图：注释单条序列 {#circular}
{{code:circular}}
{{figure:circular}}
本例显式隐藏坐标轴，让 feature 标签清楚可见。
`coord_circular()` 要求单条序列。通用 feature 可带未知链。
`geom_feature_plasmid()` 是有独立默认值的便利预设，应在需要这些默认值时显式选择。
引物和酶切图层见[识别注释](recognition.html)。

## 共线性图：对齐局部比较 {#collinear}
{{code:collinear}}
{{figure:collinear}}
`anchor` 对齐指定区间，`orientation = "auto"` 按现有证据选择支持的显示方向。
应查看方向诊断，不能直接把显示翻转解释为生物学倒位。ribbon 默认连接相邻显示行；
使用 `adjacent_only = FALSE` 或明确 `pairs` 保留其他连接。非相邻边框默认虚线，
显式样式优先。这些规则与弦图不同。[局部视图与分段](views.html)为裁剪或局部翻转
保留源坐标。
