<img class="package-logo" src="../assets/logo.webp" alt="ggchord" width="150" height="150">

面向**弦图比较、单序列圆图注释及直线共线性图**的图形语法。
ggchord 将序列骨架、区间/点连接、基因、通用 feature、标签、基因组轴和定量轨道
组合进 ggplot2 工作流。

## 阅读一张完整序列比较图 {#showcase}
{{figure:reallinear}}

四条原始序列、GFF3 注释及已有 BLAST 连接，不重编码源链方向。
[完整的数据到图形流程](real-data.html)提供筛选、报告、来源追踪及弦图视图。
另一个[完整圆图流程](plasmids.html)组合参考 feature、引物、酶切位点、堆叠与标签。

## 选择下一步 {#learn}
- [安装](installation.html)，再绘制[三种最小视图](modes.html)。
- 学习[数据契约](data.html)、[文件格式](formats.html)及[图形语法](adjustments.html)。
- 浏览[任务专题](articles.html)、[函数参考](reference.html)和[图集](gallery.html)。
- 查找当前[迁移方向](migration.html)、[FAQ](faq.html)和[版本变化](news.html)。

## 从小型可复现例子开始 {#quick-start}
{{code:setup}}
{{code:chord}}

打印 `p_chord`，再[按相同显式尺寸预览和保存](output.html)。小例子使用合成输入，
上方综合图使用原始来源记录。此文档描述基于 0.16.4 候选 API 的 0.16.5 开发版本，
正式发布版仍为 0.15.0。匹配 checkout 与旧安装版本的区别见[安装](installation.html#source)。
