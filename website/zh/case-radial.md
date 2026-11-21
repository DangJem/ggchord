这张径向图在四根辐条上分别排列四条全长序列。输入与弦图案例相同，均为 12 条 HSP 和 20 条 CDS，可以直接对照两种布局。

<a id="map"></a>

## 绘制完整图 {#map}

该图使用 12 × 10.5 英寸画布、144 dpi。脚本返回绘图对象，可传给 `ggsave()` 按相同尺寸保存。

{{example:case-radial}}

<a id="reading"></a>

## 读图 {#reading}

每条彩色骨架有自己的 bp 刻度，旁边的箭头和标签给出 CDS 位置。连接非相邻辐条时，带子向外弯曲，为绕过中间辐条留出空间。

<a id="inputs"></a>

## 来源记录与展示选择 {#inputs}

在 `layout_sequences()` 中使用随包记录的序列方向，并设置中央开口大小。在 `coord_radial_synteny()` 中用 `start_angle = 45` 旋转布局。连接带图层控制曲率和外拱，端点仍取原始坐标。

[首张径向比较图](radial-synteny.md) · [同组序列的弦图](case-chord.md)
