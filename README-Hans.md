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

### Radial · 径向辐条上的共线性

从弦图旗舰复用来源准备与精选特征，改用径向辐条骨架展示区间关系；序列身份与原生位置保持可追溯。

![从四基因组弦图派生的 Radial 旗舰](man/figures/flagradial.png)

[完整案例](website/zh/radial-synteny.md) · [绘图源码](tools/site-flagships.R)

### Sankey · 镜像双列中的来源区间

从共线性旗舰派生，保留对应的精选基因与 HSP。序列镜像排列在两列中，连接带表示已有区间关系，带宽不代表守恒流量。

![从四基因组共线性图派生的 Sankey 旗舰](man/figures/flagsankey.png)

[完整案例](website/zh/sankey.md) · [绘图源码](tools/site-flagships.R)

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

<a id="design-inspiration"></a>
## 设计灵感来源

以下项目为 ggchord 的绘图语法、序列布局与细节视图提供设计灵感，各项分别说明具体的启发方向。

### 绘图语法与序列关系表达

- [ggplot2](https://ggplot2.tidyverse.org/) — 提供分层绘图语法，以及映射、标度、图例、几何、坐标与主题的职责划分。ggchord 沿用这一模型，让用户组合独立图层，并使用自己的标度和样式。
- [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) — 为紧凑的 query/subject 字段命名提供参照，便于追踪比对记录、端点与序列标识。绘图读取已有结果，不在绘图阶段运行比对程序。
- [Circos](https://circos.ca/) · [circlize](https://github.com/jokergoo/circlize) — 为环形信息层级、克制的连接带表达、多轨道组织及关系阅读顺序提供参照，帮助组织序列骨架、注释与连接之间的层次。

### 基因箭头、注释布局与局部比较

- [clinker](https://github.com/gamcil/clinker) · [gggenes](https://github.com/wilkox/gggenes) · [genoPlotR](https://genoplotr.r-forge.r-project.org/) — 为保留链方向的基因箭头、基因簇比较和直线共线性提供清晰范例，帮助协调序列方向、特征位置与连接关系的表达。
- [gggenomes](https://thackl.github.io/gggenomes/) — 为 tidy 多表输入、focus/sync 工作流、特征堆叠、局部锚定与局部显示翻转提供参照，帮助区分数据整理和显示选择的职责。
- [DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) — 为可读的特征注释、拥挤标签处理和局部特征布局提供参照，帮助明确标签与其来源区间之间的对应关系。

### 局部放大、对齐字符与图形预览

- [ggmagnify](https://github.com/hughjonesd/ggmagnify) — 为整体图与局部放大图之间的关联表达提供参照。ggchord 围绕明确的来源区间和独立注册的序列视图发展这一表达方式。
- [ggmsa](https://bioconductor.org/packages/ggmsa/) · [ggseqalign](https://bioconductor.org/packages/ggseqalign/) — 为 DNA、蛋白对齐显示、字符着色及与 ggplot2 图形组合提供设计研究参照，尤其关注对齐列、gap 与字符细节之间的对应关系。
- [ggseqlogo](https://omarwagih.github.io/ggseqlogo/) · [Jalview](https://www.jalview.org/help/html/features/annotation.html) — 在字符视图的设计研究中，用于比较序列 logo、对齐列注释、共识与保守性的表达方式，以及逐列统计与单条序列之间的关系。
- [ggview](https://github.com/idmn/ggview) — 为按预定导出尺寸检查静态图形提供流程参照。ggchord 用 view_ggchord() 进行指定尺寸预览，用标准 ggsave() 保存。

[问题反馈](https://github.com/DangJem/ggchord/issues) · [MIT 许可证](LICENSE)
