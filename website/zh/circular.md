单序列圆图可以显示特征在一条参考序列上的位置。随包的 pUC19 记录包含这条 2,686 bp 的质粒及其注释。下面两个示例都从这份记录开始。

<a id="draw"></a>

## 画出 pUC19 骨架 {#draw}

把序列表交给 `ggchord()`，然后加上 `geom_seq()`。默认 `start_angle = 0`，起点位于十二点方向。沿顺时针方向，刻度上的坐标逐渐增大，首尾两个碱基在顶部相接。

{{example:puc-backbone}}

<a id="adapt"></a>

## 沿骨架画出九段注释区间 {#adapt}

显式传入同一记录的特征表。用 `position_feature_stack()` 分开重叠区间，用记录中的颜色或注释类型区分特征。AmpR 编码区横跨 1,626–2,486 bp。靠近序列起点的一组较短注释包括 M13 引物结合区和多克隆位点。

{{example:first-circular}}

下一步可[保留九条注释并只标出两个 CDS 的名称](annotations.md)，或[沿 pUC19 绘制 GC 含量](tracks.md)。[完整质粒图](case-circular.md)在 33 个切割位置标出 39 条限制性位点记录，[质粒图集](plasmid-catalog.md)还收录了另外十二条参考序列的示例。
