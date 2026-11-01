<img class="package-logo" src="../assets/logo.webp" alt="ggchord" width="150" height="150">

用 ggplot2 的方式绘制**弦图、单序列圆图和共线性图**。序列、连接、基因、标签和定量数据可以逐层叠加；每层的数据、配色和摆放方式都能单独设置。

## 弦图：比较多条序列 {#chord}
{{figure:realchord}}

连接表示已有比对，外侧箭头表示选出的注释。示例来自仓库中的四条原始序列。
[查看完整代码和数据筛选过程](real-data.html#chord)。

## 圆图：注释一条质粒 {#circular}
{{figure:plasmid}}

pUC19 的基因、复制起点、引物和选定酶切位点放在同一张图中。
[学习绘制质粒图](plasmids.html)。

## 共线性图：逐行查看局部关系 {#collinear}
{{figure:reallinear}}

序列按行排列，便于比较位置、方向和区间连接。
[从原始文件开始](real-data.html#linear)，或先试[三类图的简短示例](modes.html)。

## 开始使用 {#learn}
- [安装与加载](installation.html)，了解[最小输入表](data.html#tables)。
- 阅读[使用教程](articles.html)、[函数文档](reference.html)和[示例图集](gallery.html)。
- 调整[标签](labels.html)、[数据轨道](tracks.html)与[配色](styling.html)。
- 查阅[迁移指南](migration.html)、[常见问题](faq.html)和[更新日志](news.html)。

## 一个简短例子 {#quick-start}
先从随包数据中选出本例需要的列，再交给绘图函数。绘图不要求带上原始表的所有列。

{{code:setup}}
{{code:chord}}

打印 `p_chord`，再[按指定尺寸预览和保存](output.html)。当前文档对应本地 0.16.5
开发版本，正式发布版仍为 0.15.0；[安装说明](installation.html#source)解释两者区别。
