## 设计参考

ggchord 的设计参考了成熟的开源项目及科学可视化方法。这些参考涵盖绘图语法、基因组布局、注释策略和多视图组合等方面。以下分别说明各项目对 ggchord 设计过程的具体参考价值；列出参考对象不代表直接复用了其算法或实现。

### 绘图语法与 API 架构

- [ggplot2](https://ggplot2.tidyverse.org/) — 提供分层绘图语法基础，明确数据映射、几何、统计、位置、标度、图例、坐标与主题的职责。ggchord 在这一架构上扩展基因组可视化能力，支持图层组合和用户自定义样式。
- [gggenomes](https://thackl.github.io/gggenomes/) — 为基因组序列与注释的数据组织、链方向布局、特征堆叠、序列方向调整、锚点对齐和局部聚焦工作流提供重要参考。ggchord 借鉴这些领域设计理念，同时保留独立的参考坐标与图层数据契约。
- [gggenes](https://wilkox.org/gggenes/) — 为具有生物学方向的基因箭头、链方向相关的位置调整，以及基因组区间长度与图形物理尺寸的区分提供参考。

### 环形与比较基因组可视化

- [Circos](https://circos.ca/) — 为环形基因组信息的视觉组织提供灵感，包括序列扇区、多层注释和表示基因组区间关系的连接带。
- [circlize](https://github.com/jokergoo/circlize) — 为灵活的环形轨道组织、环层排列和多层基因组可视化提供设计参考。
- [clinker](https://github.com/gamcil/clinker) — 为基因簇比较、方向性基因箭头，以及基因组相似性和共线性关系的清晰表达提供视觉参考。
- [genoPlotR](https://genoplotr.r-forge.r-project.org/) — 为线性比较基因组布局、序列区域对齐、基因组坐标表达和多个基因组片段之间的关系展示提供参考。

### 基因组注释与多尺度可视化

- [DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) — 为基因组特征展示、注释几何、标签位置以及高密度注释区域中的碰撞处理提供参考。
- [ggmagnify](https://github.com/hughjonesd/ggmagnify) — 为整体图与局部放大视图之间的关联表达提供参考。ggchord 进一步围绕明确的来源区间和独立构建的子视图实现跨视图连接。
- [ggmsa](https://bioconductor.org/packages/ggmsa/) — 为核酸与氨基酸比对可视化、紧凑残基单元格及化学性质配色提供参考。部分残基颜色值经过适配并在 [NOTICE](../inst/NOTICE) 中注明来源；序列注册与图形渲染由 ggchord 独立实现。
- [ggview](https://github.com/idmn/ggview) — 为按照预定输出尺寸预览科研图形的工作流提供参考。ggchord 通过 `view_ggchord()` 提供尺寸感知的预览，并继续使用标准 `ggsave()` 保存图形。

### 数据规范与相关研究

- [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) — 其成熟的 query/subject 术语为 ggchord 的比对数据字段规范提供参考。ggchord 负责可视化已有比对结果，不在绘图过程中执行序列比对。
- [trackViewer](https://bioconductor.org/packages/trackViewer/) — 作为基因组信号轨道、区间测量及定量注释展示的研究参考对象。
- [ggseqalign](https://bioconductor.org/packages/ggseqalign/)、[ggseqlogo](https://omarwagih.github.io/ggseqlogo/)、[Jalview](https://www.jalview.org/) 和 [ESPript](https://espript.ibcp.fr/ESPript/ESPript/) — 作为序列比对可视化、序列 Logo、残基注释及发表级序列图形的相关研究与比较对象。列入此处不代表已完成工具集成或采用其图形渲染算法。
