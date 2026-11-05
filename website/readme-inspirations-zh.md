## 设计参考与致谢

ggchord 的图形语法、基因组表达和预览流程参考了以下开源项目与公共工具：

- [ggplot2](https://ggplot2.tidyverse.org/)：以图层组合图形，区分数据映射、标度、几何、坐标与主题。
- [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/)：借鉴 query/subject 表格字段命名，便于追踪比对端点与序列标识。
- [Circos](https://circos.ca/) 和 [circlize](https://github.com/jokergoo/circlize)：参考环形信息层级、连接带表达与多轨道的阅读秩序。
- [clinker](https://github.com/gamcil/clinker)、[gggenes](https://github.com/wilkox/gggenes) 和 [genoPlotR](https://genoplotr.r-forge.r-project.org/)：参考基因箭头、基因簇与线性共线性比较的可读表达。
- [DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) 和 [gggenomes](https://thackl.github.io/gggenomes/)：参考注释布局、多表输入、特征堆叠及局部比较工作流。
- [ggview](https://github.com/idmn/ggview)：参考按实际导出尺寸预览的流程；ggchord 用 `view_ggchord()` 预览，用标准 `ggsave()` 保存。

这些项目为设计与比较提供了重要参考。具体可用功能、参数及运行依赖以 ggchord 的函数文档和 `DESCRIPTION` 为准。
