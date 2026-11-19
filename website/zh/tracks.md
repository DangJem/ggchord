<a id="shared-band"></a>

## 在同一高度轴叠加柱与折线 {#shared-band}

随包 bedGraph 记录 pUC19 各窗口的 GC 百分数。柱横跨原窗口，折线连接窗口中点。用 `track_band("gc")` 声明一次物理轨带，两种几何都用 `track = "gc"` 选取它。两层高度为 60 的位置都代表 60% GC。

{{example:gc-shared-band}}

<a id="domain"></a>

## 分别控制物理宽度、数值和颜色 {#domain}

轨带的 `width` 决定它在纸上占多宽。`scale_track_value_continuous()` 的 `limits` 决定数值轴范围。更改带宽不影响窗口的 GC 百分数。若想用颜色表示区间数值，可在 tile 层映射 `aes(fill = value)`。有不同单位的第二种度量时，另设一条带及其高度轴。

若需要从区间计算覆盖率再作分组比较，见[逐窗统计与分面](tracks-statistics.md)。
