从随包 pUC19 bedGraph 已经记录的 GC 百分比开始。定量轨道在序列旁为这些数值提供一条高度轴。下面使用可选的 dplyr 包整理表格。

<a id="domain"></a>

## 每个实测窗口画一根柱 {#domain}

用 `read_bedgraph()` 读取已安装的文件。前三个 100 bp 窗口的 GC 含量分别为 60%、57% 和 59%。每个窗口的值画一根柱。图层自动创建轨带，用 `layout_track()` 设置厚度，用 `scale_y_continuous()` 把百分比轴设为 0–100。第一根柱从第 1 bp 延伸至第 100 bp。

{{example:gc-bars}}

<a id="shared-band"></a>

## 沿同一组窗口值叠加折线 {#shared-band}

计算各窗口的中点。保留柱图层，再加一层折线，以中点为横坐标、实测值为纵坐标。用 `add_track()` 将两层放进同一轨带，共用 Y 标度。折线会经过每根柱顶边的中点。

{{example:gc-shared-band}}

轨带的 `thickness` 预留 15 mm 的绘制空间，`scale_y_continuous()` 把测量值映射为这段空间中的高度。加厚轨带可让一个百分点的差异更易辨认，测量值仍保持原样。下一步可[从比对区间计算覆盖率](tracks-statistics.md)，或[查看 bedGraph 的原始字段](import-signals.md)。
