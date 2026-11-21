这是首页圆图预览所链接的完整 pUC19 图。它将来源注释与酶切位点共同放在 2,686 bp 的双线骨架上。

<a id="map"></a>

## 绘制完整图 {#map}

该图使用 10 × 8 英寸画布、110 dpi。脚本返回绘图对象，可传给 `ggplot2::ggsave()` 按相同尺寸保存。

{{example:case-circular}}

<a id="reading"></a>

## 读图 {#reading}

九条特征均保留参考坐标和链。沿骨架刻度定位特征，再从圆环外侧读取酶切标签。47 条位点记录对应 32 个位置，同一位置的酶名并入一个标签。这份记录没有独立参考引物。

<a id="inputs"></a>

## 来源记录与展示选择 {#inputs}

`flagship_circular_example$references[["pUC19"]]` 中保存完整长度和三张用于绘图的表。特征层设置分道位置，标签层按图层标识跟随所用实体；代码会在参考记录包含引物行时加入引物层。

[只绘制九条特征](circular.md) · [对照十三张完整参考图](plasmid-catalog.md)
