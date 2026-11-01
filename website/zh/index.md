<div class="home-intro">
<img class="package-logo" src="../assets/logo.webp" alt="ggchord" width="150" height="150">
<p class="eyebrow">GENOMIC VISUALIZATION · R / GGPLOT2</p>
<p class="home-lead">让序列关系，清晰可见。</p>
<p>用 ggplot2 图层组合<strong>弦图、单序列圆图与共线性图</strong>。从已有的序列、注释和比对出发，逐步构建可解释、可复现的基因组图形。</p>
<div class="home-actions"><a class="primary-action" href="installation.html">绘制第一张图 →</a><a href="articles.html">查找任务教程</a><a href="reference/index.html">查阅函数</a></div>
</div>

## 先确定你想回答的问题 {#choose}
| 你想了解什么？ | 选择 | 从哪里开始 |
| --- | --- | --- |
| 多条序列之间有哪些比对关系？ | 弦图，强调总体连接与分布 | [真实序列比较](real-data.html) |
| 一条圆形序列上有哪些注释？ | 圆图，保留完整序列与方向 | [pUC19 注释流程](plasmids.html) |
| 某段基因在不同序列上的位置与方向如何？ | 局部共线性图，让对应区间靠近阅读 | [比较真实 CDS](homology.html) |

不是每个问题都适合圆形布局。[选图指南](modes.html)解释输入需求、阅读方式及适用边界。

## 弦图 · 看总体关系 {#chord}
{{figure:realchord}}

四条噬菌体序列的已有比对。序列用类别色区分，连接按真实 `pident` 连续配色；外侧箭头只显示明确筛选的 CDS。先看连到哪里，再读相似度，交叉并不直接证明重排。
[读取原始 FASTA、GFF3 和 BLAST，重现这张图 →](real-data.html#read)

## 圆图 · 读单序列注释 {#circular}
{{figure:plasmid}}

pUC19 的 CDS、复制起点与启动子使用不同类别色；引物和酶切位点使用独立线与标签。图上有意只显示三种酶，完整匹配结果仍保留。
[从随包 DNA 到完整质粒图 →](plasmids.html)

## 共线性图 · 比较局部位置 {#collinear}
{{figure:homologs}}

pUC19 与 pETDuet-1 中编码 DNA 完全相同的 AmpR CDS。两侧各保留 200 bp，坐标轴显示源位置，同色箭头表示有序列证据的对应 CDS。连带表示已声明的对应范围，未重新计算比对。
[核对 CDS、建立关系并聚焦邻域 →](homology.html)

## 按你的起点继续 {#learn}
<div class="route-grid">
<a class="route-card" href="installation.html"><strong>第一次使用</strong><span>安装 → 最小输入 → 首图 → 读图 → 保存</span></a>
<a class="route-card" href="formats.html"><strong>已有自己的数据</strong><span>选择读取器，确认标识、长度和坐标单位</span></a>
<a class="route-card" href="styling.html"><strong>改进已有图形</strong><span>配色、图例、标签、局部视图与最终尺寸</span></a>
<a class="route-card" href="reference/index.html"><strong>查找精确参数</strong><span>完整函数签名、参数、返回值与示例</span></a>
</div>

## 用随包数据开始 {#quick-start}
先绘制四条随包序列的骨架，检查标识和完整长度：

```r
library(ggchord)
data(seq_data_example)
sequences <- seq_data_example[c("accver", "length")]
ggchord(seq_data = sequences) +
  geom_seq() +
  geom_seq_label() +
  coord_chord()
```

[开始使用](installation.html#input)给出独立可运行的准备代码，只将序列两列、基因四列和连接六列传入绘图函数。无需先学习全部图层或整理历史元数据。

已有旧版本代码？先看[迁移指南](migration.html)。当前文档对应本地 0.16.5 开发版；正式发布版仍为 0.15.0，安装差异见[版本说明](installation.html#source)。
