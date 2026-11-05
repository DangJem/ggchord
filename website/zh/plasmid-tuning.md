对照默认圆图，分别改骨架样式和整图旋转。

## 默认双线骨架 {#baseline}

{{code:firstcircle}}

{{figure:firstcircle}}

## 切换为单线骨架 {#style}

{{code:frontcirclestyle}}

{{figure:frontcirclestyle}}

`seq_style = "single"` 把双线改为单线。中心名称和 bp 刻度仍可用于定位。

## 转动圆图 {#rotation}

{{code:firstfeatures,frontcirclerotation}}

{{figure:frontcirclerotation}}

默认 `rotation = 90` 把序列起点放在 12 点钟，坐标顺时针增加。这里设置 `rotation = 180`，起点移到 9 点钟；圆周刻度和复制起点一起旋转。原始注释区间仍是同一段 bp。下一节[整理自己的序列和注释](plasmid-data.html)。
