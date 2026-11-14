<!-- Generated from website/readme-template.md and website/readme-text.json. -->
<p align="center">
  <img src="man/figures/logo.png" alt="ggchord" width="280">
</p>

[English](README.md) · [现代汉语](README-Hans.md)

# ggchord：序列、注释与关系

基于 ggplot2 的 R 可视化包：比较多序列关系、组织单序列注释、展示 HSP 比对，并组合全局与局部视图。

[网站](https://dangjem.github.io/ggchord/) · [使用说明](website/zh/modes.md) · [示例图集](website/zh/gallery.md) · [问题反馈](https://github.com/DangJem/ggchord/issues)

<a id="modules"></a>
## 四个绘图模块

| 模块 | 主要用途 | 视图 |
| --- | --- | --- |
| [I · 多序列关系](website/zh/multi-sequence.md) | 比较序列区间、基因顺序与已有连接关系。 | Chord · Collinear · Radial · Sankey |
| [II · 单序列圆图](website/zh/circular.md) | 在一条序列上组织特征、引物和限制位点。 | Circular |
| [III · HSP 点图](website/zh/dotplot.md) | 在 query 与 subject 双轴上查看比对区间。 | Dotplot |
| [IV · 多视图组合](website/zh/fusion.md) | 用明确的来源区间连接独立图形。 | Fusion |

多序列模块包含四种 Coord；Circular、Dotplot 各有专属坐标系，Fusion 负责组合独立图形。图层、映射、标度、局部位置与主题各司其职，沿用 ggplot2 的绘图语法。

<a id="examples"></a>
## 真实数据综合示例

### Chord · 四条基因组的关系总览

混合曲率骨架、保留真实链方向的 CDS 箭头、来源标签与 Identity 连接带组合展示。特征与 HSP 按显示需要从原始记录中选取。

![四条基因组的弦图旗舰：基因标签与比对连接带](man/figures/showcasechord.png)

[完整案例](website/zh/chord.md) · [四种关系视图转换](website/zh/multi-sequence-views.md) · [绘图源码](tools/documentation-examples.R)

### Collinear · 顺序、方向与跨行连接

在共同的 bp 标尺下比较四条基因组，骨架两侧展示精选 CDS。跨越非相邻显示行的连接带使用虚线轮廓。

![四条基因组的共线性旗舰：分链基因与跨行连接带](man/figures/showcasecollinear.png)

[完整案例](website/zh/collinear.md) · [对齐与排列](website/zh/collinear-layout.md) · [绘图源码](tools/documentation-examples.R)

### Circular · pETDuet-1 密集注释

完整圆图包含 20 个特征、5 个参考引物，以及声明酶集合中的 58 个合并切点位置。序列位置与注释身份均来自原始输入。

![pETDuet-1 圆图旗舰：完整特征、引物与限制位点](man/figures/showcaseplasmid.png)

[完整案例](website/zh/circular-dense-case.md) · [13 个完整参考图](website/zh/plasmid-catalog.md) · [绘图源码](tools/documentation-examples.R)

### HSP 矩阵与多视图组合

<table>
  <tr>
    <th>III · HSP 点图</th>
    <th>IV · 多视图组合</th>
  </tr>
  <tr>
    <td width="50%"><a href="man/figures/flagdotplot.png"><img src="man/figures/flagdotplot.png" alt="四序列、32 条原始 HSP 的双轴矩阵" width="460"></a></td>
    <td width="50%"><a href="man/figures/flagfusion.png"><img src="man/figures/flagfusion.png" alt="连接 Chord、Collinear、Radial 与局部放大的组合图" width="460"></a></td>
  </tr>
  <tr>
    <td>32 条原始 HSP；颜色表示 Identity，线宽表示原记录中的比对长度。</td>
    <td>在同一组合中关联全局关系、共线性比较与来源区间的局部放大。</td>
  </tr>
  <tr>
    <td><a href="website/zh/dotplot-start.md">双轴、方向与映射 →</a></td>
    <td><a href="website/zh/fusion-details.md">局部放大与字符细节 →</a></td>
  </tr>
</table>

图像的输入、显示选择与绘图定义均可追溯：
[输入文件](examples/) · [Dotplot / Fusion 源码](tools/site-flagships.R) · [数据来源与归属](website/ASSET-SOURCES.md)

<a id="installation"></a>
## 安装

R ≥ 4.1.0 · ggplot2 ≥ 4.0.0

```r
install.packages("remotes")
remotes::install_github("DangJem/ggchord", dependencies = NA)
```

下面的最小示例只需加载 ggchord。完整数据整理示例会在代码旁声明 dplyr、purrr、tidyr 等可选依赖。

<a id="quick-start"></a>
## 绘制第一张图

使用随包表格绘制序列骨架、比对连接带和序列名称，无需下载文件或运行比对程序。

```r
library(ggchord)

p <- ggchord(seq_data_example, ribbon_data_example) +
  geom_seq() +
  geom_link_ribbon() +
  geom_seq_label() +
  coord_chord()

p
```

这段最小示例展示图层组合方式。上方综合图进一步加入了对应的来源注释、显示选择和样式。

<a id="output"></a>
## 预览与保存

按相同的明确宽高预览并保存，让标签按预定成图尺寸排布。

```r
view_ggchord(p, width = 10, height = 8, units = "in")
ggplot2::ggsave("ggchord.png", p, width = 10, height = 8, units = "in", dpi = 144)
```

<a id="documentation"></a>
## 说明与资源

- [绘图基本概念](website/zh/modes.md) · [输入数据表](website/zh/data.md) · [文件导入](website/zh/formats.md)
- [标签布局](website/zh/labels.md) · [外观与主题](website/zh/styling.md) · [图形输出](website/zh/output.md)
- [诊断与排错](website/zh/diagnostics.md) · [函数参考源码](man/) · [变更记录](NEWS.md)

<details>
<summary>设计参考与致谢</summary>

ggchord 的设计与比较参考了以下项目和公共工具：

- [ggplot2](https://ggplot2.tidyverse.org/) — 分层绘图语法。
- [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) — query/subject 字段命名。
- [Circos](https://circos.ca/) · [circlize](https://github.com/jokergoo/circlize) — 环形信息层级与轨道表达。
- [clinker](https://github.com/gamcil/clinker) · [gggenes](https://github.com/wilkox/gggenes) · [genoPlotR](https://genoplotr.r-forge.r-project.org/) — 基因箭头与共线性比较。
- [DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) · [gggenomes](https://thackl.github.io/gggenomes/) — 注释布局与局部比较。
- [ggview](https://github.com/idmn/ggview) — 按导出尺寸预览图形。

</details>

[问题反馈](https://github.com/DangJem/ggchord/issues) · [MIT 许可证](LICENSE)
