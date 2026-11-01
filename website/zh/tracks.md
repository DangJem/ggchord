## 从实际 DNA 计算数据轨道 {#domain}
先运行[共享准备代码](index.html#quick-start)。本例读取随包 pUC19，每 100 bp 计算
一次 GC 比例；末尾不足 100 bp 的窗口按实际碱基数计算。

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
先运行[圆图例子](modes.html#circular)得到实际注释 `features`，再继续：

{{code:multitrack}}
{{figure:multitrack}}

柱和线都来自同一组实际 GC 值，线使用三个窗口的滑动均值；两端没有完整窗口时省略，
不填造观测值。折线位置是整数窗口中点。注释类型带有独立 ID 和配色。
轨道放在圆外并按顺序排列；小半径下过宽的内侧带可能折叠，此时会明确报错。
