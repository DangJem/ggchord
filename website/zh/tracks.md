<div class="task-intro"><p>从真实 DNA 计算 GC 比例，再逐步叠加柱、滑动均值和注释带。读者应能说清每个轨道的单位、窗口与基线。</p></div>

## 本页准备 {#prepare}
在新的 R 会话中先运行下方代码，再按页面顺序继续。配色和排版是可选的显示设置，输入只选择需要的列。

{{prepare:docstyle,plasmiddata}}

## 从实际 DNA 计算数据轨道 {#domain}
本例读取随包 pUC19，每 100 bp 计算
一次 GC 比例；末尾不足 100 bp 的窗口按实际碱基数计算。

{{code:trackdata}}

`bins` 的每一行是一段实际序列窗口；先检查起止位置及 `gc`，再把数值映射到轨道高度与颜色。

### 绘制 GC 比例

{{code:tracks}}
{{figure:tracks}}

值域为 0–1，单位是比例。换成其他测量值时，要声明对应单位和有效范围。
`position_track()` 决定放置与物理宽度，`scale_track_value_continuous()` 决定数值映射，
二者不能互相替代。

## 同一 ID 共享同一条轨道 {#sharing}
同名轨道上的各层共享数值域，bar 的基线也参与范围计算；它们必须声明一致的宽度和
放置方式。颜色标度按轨道 ID 独立设置。tile 表示区间类别，line 连接观测位置，
bar 表示相对基线的区间值。不要把相同数值域误认为必须位于同一条轨道。

## 检查裁剪和数值范围 {#diagnostics}
`export_ggchord_layout()` 的 `track_registry` 与 `track_diagnostics` 记录值域、
物理尺寸、原生位移及裁剪情况。被裁剪或省略的行不等于测量值为零。
局部视图使用 `focus_ggchord_track()` 并共享相同窗口。

## 叠加柱、滑动均值和注释带 {#multiple}
本页准备代码已创建实际注释 `features`，完成上面的 GC 计算后继续：

{{code:multitrack}}
{{figure:multitrack}}

柱和线都来自同一组实际 GC 值，线使用三个窗口的滑动均值；两端没有完整窗口时省略，
不填造观测值。折线位置是整数窗口中点。注释类型带有独立 ID 和配色。
轨道放在圆外并按顺序排列；小半径下过宽的内侧带可能折叠，此时会明确报错。

## 完成后检查 {#check-result}
bins 的 gc 在 0–1 之间，每行对应一个实际窗口；最后窗口可短于 100 bp。连续色表达数值，类别带只表达注释类型，不能用它比较 GC 大小。

<div class="next-step">继续：<a href="output.html">预览和保存图形 →</a></div>

## 保存本页结果 {#save-result}
在工作目录中保存刚刚创建的图；预览时使用相同宽高。

```r
ggsave("pUC19-tracks.png", p_multitrack, width = 9, height = 8, units = "in", dpi = 144)
```
