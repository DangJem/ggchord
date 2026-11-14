## 从实际 DNA 计算数据轨道 {#domain}
本例读取随包 pUC19，每 100 bp 计算
一次 GC 比例；末尾不足 100 bp 的窗口按实际碱基数计算。

{{code:trackdata,frontgctable}}

{{output:frontgctable}}

`bins` 的每一行是一段实际序列窗口；先检查起止位置及 `gc`，再把数值映射到轨道高度与颜色。

### 绘制 GC 比例 {#draw-gc}

{{code:tracks}}
{{figure:tracks}}

值域为 0–1，单位是比例。换成其他测量值时，要声明对应单位和有效范围。
`position_track()` 决定放置与物理宽度，`scale_track_value_continuous()` 决定数值映射，
二者不能互相替代。

## 同一 ID 共享同一条轨道 {#sharing}
同名轨道上的各层共享数值域，bar 的基线也参与范围计算；它们必须声明一致的宽度和
放置方式。颜色标度按轨道 ID 独立设置。tile 表示区间类别，line 连接观测位置，
`geom_track_col()` 表示给定区间值；`geom_track_bar()` 计数同位置观测，基线仅允许 0。bin 统计才要求显式 binwidth/breaks。不要把相同数值域误认为必须位于同一条轨道。

## 检查裁剪和数值范围 {#diagnostics}
`export_ggchord_layout()` 的 `track_registry` 与 `track_diagnostics` 记录值域、
物理尺寸、原生位移及裁剪情况。被裁剪或省略的行不等于测量值为零。
局部视图使用 `focus_ggchord_track()` 并共享相同窗口。

## 叠加柱与滑动均值 {#multiple}
沿用上面的 GC 窗口，加一条三窗口滑动均值线。

{{code:multitrack}}
{{figure:multitrack}}

柱和线都来自同一组实际 GC 值，线使用三个窗口的滑动均值；两端没有完整窗口时省略，
不填造观测值。折线位置是整数窗口中点；柱与线共用同一个轨道 ID。

## 用颜色绘制 tile 轨道 {#tile}

沿用本页从真实 DNA 得到的 `track_seq` 与 `bins`。tile 要求 `accver, start, end`，并不要求数值高度；这里把同一个 `gc` 映射到 `track_fill`。

{{code:tracktile}}
{{figure:tracktile}}

与柱图使用相同 9 × 8 英寸尺寸和 100 bp 窗口。柱高表示数值；tile 的固定带宽不随 GC 改变，颜色才表示 0–1 的 GC 比例。它不是另一组实验观测，也不会自动平均重叠区间。此图比较编码方式，**不是单参数几何对照**。

| 参数或映射 | 默认值 | 本例选择与限制 |
| --- | --- | --- |
| `position` | 必须显式提供 | `position_track("gc_colour")` 声明独立轨道 |
| `width` | `position_track()` 为 4 mm | 本例 10 mm；是物理带宽 |
| `side` / `order` | outside / 1 | 多轨道按 order 排列，越小越靠近骨架 |
| `gap` | 1 mm | 到前一带的物理间隔，不是序列缺口 |
| `track_fill` | 不自动推断 | 在 `aes()` 映射 gc；固定颜色写 `fill` |
| `bounds` | `"error"` | 越界默认报错；clip 会报告裁剪 |
| `alpha` / `linewidth` | 1 / 0.25 | 固定透明度和边框线宽；默认边框 colour 为 NA |

连续颜色用 `scale_track_fill_gradientn()`；类别用 `scale_track_fill_manual()`。二者都要填写对应轨道 ID。不要把比例误标为百分数。若要叠加本页柱图，把该 tile 图层加到 `p_tracks`，保留不同 ID：这是两条独立轨道。相同 ID 必须使用一致的 side/order/width/gap。

试把 tile 的 width 从 10 改为 4 mm，其他设置不动；数值颜色不变，只有物理带宽变化。自备数据必须声明窗口、单位、NA 和重叠区间的处理。源坐标越界先查输入，不能用 clip 掩盖错误。更多见 [`geom_track_tile()`](reference/geom_track_tile.html) 和 [`position_track()`](reference/position_track.html)。

## 只改变 tile 的物理宽度 {#tile-width}

回看上面的 10 mm 色带，保持 9 × 8 英寸、输入窗口、颜色和数值域完全相同，仅把 `position_track(width)` 改为 4 mm：

{{code:comparetrackthin}}
{{figure:comparetrackthin}}

颜色仍表示同一批 GC 比例，变窄的是物理带宽。`geom_track_tile()` 负责区间色块，Position 负责带宽；它们不改变窗口里的碱基组成。把 4 改为 6 mm 操作示例一次，保留同一导出尺寸。柱图的数值高度另由 `track_value` 和 value scale 控制，不能用缩窄带宽代替数值归一化。

## 检查要点 {#check-result}
bins 的 gc 在 0–1 之间，每行对应一个实际窗口；最后窗口可短于 100 bp。柱高与折线共用 0–1 的比例尺度。

## 保存本页结果 {#save-result}
在工作目录中保存刚刚创建的图；预览时使用相同宽高。

{{code:courseblockv}}
