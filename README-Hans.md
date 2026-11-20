<!-- Generated from website/readme-template.md and website/readme-text.json. -->
<p align="center">
  <img src="man/figures/logo.png" alt="ggchord" width="280">
</p>

[English](README.md) · [现代汉语](README-Hans.md)

# ggchord：序列、注释与关系

基于 ggplot2 的 R 可视化包：比较多序列关系、组织单序列注释、展示 HSP 比对，并组合全局与局部视图。

> [!IMPORTANT]
> 此 README 对应开发中的 v0.50.0 源码。请从当前检出目录安装后运行示例。

[网站](https://dangjem.github.io/ggchord/) · [教程](vignettes/ggchord_guide_hans.Rmd) · [问题反馈](https://github.com/DangJem/ggchord/issues)

<a id="modules"></a>
## 四个绘图模块

| 模块 | 主要用途 | 视图 |
| --- | --- | --- |
| I · 多序列关系 | 比较序列区间、基因顺序与已有连接关系。 | Chord · Collinear · Radial · Sankey |
| II · 单序列圆图 | 在一条序列上组织特征、引物和限制位点。 | Circular |
| III · HSP 点图 | 在 query 与 subject 双轴上查看比对区间。 | Dotplot |
| IV · 多视图组合 | 用明确的来源区间连接独立图形。 | Fusion |

多序列模块包含四种 Coord；Circular、Dotplot 各有专属坐标系，Fusion 负责组合独立图形。图层、映射、标度、局部位置与主题各司其职，沿用 ggplot2 的绘图语法。

<a id="examples"></a>
## 真实数据综合示例

### I · 多序列关系

<table>
  <tr>
    <td width="50%"><strong>Chord · 四条基因组的关系总览</strong><br><a href="man/figures/showcasechord.png"><img src="man/figures/showcasechord.png" alt="四条基因组的弦图旗舰：基因标签与比对连接带" width="460"></a><br>比较四条基因组参考。基因箭头显示链向，连带按原始比对一致性着色。<br><a href="examples/make_readme_plots.R">完整案例 →</a></td>
    <td width="50%"><strong>Collinear · 顺序、方向与跨行连接</strong><br><a href="man/figures/showcasecollinear.png"><img src="man/figures/showcasecollinear.png" alt="四条基因组的共线性旗舰：分链基因与跨行连接带" width="460"></a><br>同一组源区间使用共同的碱基位置比例，按显式声明的显示方向排列。<br><a href="examples/make_readme_plots.R">完整案例 →</a></td>
  </tr>
</table>

<table>
  <tr>
    <td width="50%"><strong>Radial · 径向辐条上的共线性</strong><br><a href="man/figures/flagradial.png"><img src="man/figures/flagradial.png" alt="从四基因组弦图派生的 Radial 旗舰" width="460"></a><br>同一组观察投影到径向辐条，保留参考身份与源区间端点。<br><a href="examples/make_readme_plots.R">完整案例 →</a></td>
    <td width="50%"><strong>Sankey · 镜像双列中的来源区间</strong><br><a href="man/figures/flagsankey.png"><img src="man/figures/flagsankey.png" alt="从四基因组共线性图派生的 Sankey 旗舰" width="460"></a><br>各参考在查询列和目标列各显示一次，连带颜色沿用其他比较视图的一致性标度。<br><a href="examples/make_readme_plots.R">完整案例 →</a></td>
  </tr>
</table>

### II · 单序列圆图

<a href="man/figures/showcaseplasmid.png"><img src="man/figures/showcaseplasmid.png" alt="pETDuet-1 圆形注释图，含区间特征和酶切识别位置" width="460"></a>

pETDuet-1 图组合随包区间注释、选定的标签和单切酶的识别位置。

[完整案例](examples/make_readme_plots.R)

### III · HSP 点图

<a href="man/figures/flagdotplot.png"><img src="man/figures/flagdotplot.png" alt="四序列、32 条原始 HSP 的双轴矩阵" width="460"></a>

显示 32 条原始 HSP，颜色表示原始一致性，线宽表示原始比对长度。

[完整案例](examples/make_readme_plots.R)

### IV · 多视图组合

<a href="man/figures/flagfusion.png"><img src="man/figures/flagfusion.png" alt="由两个面板组成的源区间全局图与局部放大图" width="460"></a>

弦图总览与独立的共线性细节图以不同尺度展示同一来源区间。

[完整案例](examples/make_readme_plots.R)

图像的输入、显示选择与绘图定义均可追溯：
[输入文件](examples/) · [绘图源码](examples/make_readme_plots.R) · [数据制备](data-raw/) · [数据来源与归属](website/ASSET-SOURCES.md)

<a id="get-started"></a>
## 开始使用

<a id="installation"></a>
### 安装

R ≥ 4.1.0 · ggplot2 ≥ 4.0.0

从 [CRAN](https://CRAN.R-project.org/package=ggchord) 安装稳定版：

```r
install.packages("ggchord")
```

在当前检出目录的根目录安装开发版源码。

```r
install.packages("remotes")
remotes::install_local(".", dependencies = NA)
```

入门图使用随包数据和 ggplot2。旗舰图复现脚本还使用 dplyr。

<a id="quick-start"></a>
### 绘制第一张图

使用随包表格绘制序列骨架、比对连接带和序列名称，无需下载文件或运行比对程序。

```r
library(ggchord)
library(ggplot2)

p <- ggchord(seq_data_example) +
  geom_seq() +
  geom_link_ribbon(aes(fill = pident), data = ribbon_data_example) +
  geom_seq_label() +
  coord_chord()

p
```

观察层各自接收数据。切换坐标系时，序列参考定义保持一致。

<a id="output"></a>
### 预览与保存

按相同的明确宽高预览并保存，让标签按预定成图尺寸排布。

```r
view_ggchord(p, width = 10, height = 8, units = "in")
ggplot2::ggsave("ggchord.png", p, width = 10, height = 8, units = "in", dpi = 144)
```

<a id="design-inspiration"></a>
## 设计灵感来源

以下条目区分已记录并落实到设计的参考，以及仅用于能力调研的项目。具体依据见[设计来源核对](tools/site-audit/design-inspiration.md)。

### 绘图语法与序列关系表达

- [ggplot2](https://ggplot2.tidyverse.org/) — 提供分层绘图语法，以及映射、标度、图例、几何、坐标与主题的职责划分。ggchord 沿用这一模型，让用户组合独立图层，并使用自己的标度和样式。
- [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) — 为紧凑的 query/subject 字段命名提供参照，便于追踪比对记录、端点与序列标识。绘图读取已有结果，不在绘图阶段运行比对程序。
- [Circos](https://circos.ca/) · [circlize](https://github.com/jokergoo/circlize) — Circos 为环形信息层级和连接带表达提供设计参照；circlize 用于多轨道与显式 track/ring 契约的能力比较。

### 基因箭头、注释布局与局部比较

- [clinker](https://github.com/gamcil/clinker) · [gggenes](https://github.com/wilkox/gggenes) · [genoPlotR](https://genoplotr.r-forge.r-project.org/) — 为保留链方向的基因箭头、基因簇比较和直线共线性提供清晰范例，帮助协调序列方向、特征位置与连接关系的表达。
- [gggenomes](https://thackl.github.io/gggenomes/) — 为 tidy 多表输入、focus/sync 工作流、特征堆叠、局部锚定与局部显示翻转提供参照，帮助区分数据整理和显示选择的职责。
- [DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) — 为可读的特征注释、拥挤标签处理和局部特征布局提供参照，帮助明确标签与其来源区间之间的对应关系。

### 局部放大、对齐字符与图形预览

- [ggmagnify](https://github.com/hughjonesd/ggmagnify) — 为整体图与局部放大图之间的关联表达提供参照。ggchord 围绕明确的来源区间和独立注册的序列视图发展这一表达方式。
- [ggmsa](https://bioconductor.org/packages/ggmsa/) — v0.20.0 将其列为字符视图调研对象；v0.28.0 进一步实际核对其残基格和化学性质配色，并用于 geom_residue() 的紧凑背景与默认颜色标度。保留独立字形绘制、真实坐标注册和用户标度优先规则。
[ggseqalign](https://www.bioconductor.org/packages/release/bioc/html/ggseqalign.html)、[ggseqlogo](https://omarwagih.github.io/ggseqlogo/)、[Jalview](https://www.jalview.org/help/html/features/annotation.html) 和 [ESPript](https://espript.ibcp.fr/ESPript/ESPript/) 仅列为对齐字符、序列 logo、列注释或静态输出的调研比较对象；当前未实现这些工具的集成、logo、共识或保守性图层。历史记录不足以证明它们分别影响了现有布局算法。
- [ggview](https://github.com/idmn/ggview) — 为按预定导出尺寸检查静态图形提供流程参照。ggchord 用 view_ggchord() 进行指定尺寸预览，用标准 ggsave() 保存。

[问题反馈](https://github.com/DangJem/ggchord/issues) · [MIT 许可证](LICENSE)
