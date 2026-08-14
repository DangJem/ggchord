# ggchord：基于 ggplot2 的多序列弦图可视化工具

🌐
语言切换：【[现代汉语（Hans）](https://dangjem.github.io/ggchord/README-Hans.md)
\| [English](https://dangjem.github.io/ggchord/README.md)】

## 概述

`ggchord` 是一个基于 `ggplot2` 的 R
语言包，采用**分层的图形语法**将多序列数据绘制为直观的弦图。与“一个大函数”不同，你通过叠加图层来构建图形：[`ggchord()`](https://dangjem.github.io/ggchord/reference/ggchord.md)
负责提供数据与全局选项，每个 `geom_*`
图层负责绘制一类元素（序列弧线、比对连接带、基因注释、坐标轴、标签）。每个图层都有合理的默认值，因此**一行代码即可画出完整的弦图**，需要精细控制时也可以分别微调每个图层。

该包是通用的多序列比较工具，可用于序列比较、基因邻域分析、噬菌体-宿主关系、泛基因组区块、共线性分析等——你只需要准备三张规整的数据表。

## 主要功能

- **真正的 `ggplot2` 分层风格**：数据与全局选项交给
  [`ggchord()`](https://dangjem.github.io/ggchord/reference/ggchord.md)，每个
  `geom_*` 图层管理自己的布局参数。
- **开箱即用**：`ggchord(data) + geom_seq() + geom_ribbon() + geom_gene() + geom_axis()`
  即可得到完整图形。
- **多序列支持**：可同时展示两条、三条、四条或更多序列。
- **灵活的参数**：支持单值、向量、命名向量与列表，可按序列或链方向分别设置。
- **与 ggplot2
  生态无缝衔接**：`theme()`、`scale_*()`、`ggsave()`、`ggplot_build()`、[`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html)
  均可使用。

## 安装

`ggchord` 需要 R（≥ 4.1.0）和 `ggplot2`（≥ 4.0.0）。

``` r

install.packages("ggplot2")   # 如需要
```

从 CRAN 安装：

``` r

install.packages("ggchord")
```

从 GitHub 安装（开发版）：

``` r

devtools::install_github("DangJem/ggchord")
```

## 快速开始

包内自带三份示例数据，直接运行以下代码即可：

``` r

library(ggchord)

# 加载内置示例数据
data(seq_data_example)
data(ribbon_data_example)
data(gene_data_example)

# 像 ggplot2 一样叠加图层
ggchord(
  seq_data     = seq_data_example,
  ribbon_data  = ribbon_data_example,
  gene_data    = gene_data_example
) +
  geom_seq() +      # 序列弧线
  geom_ribbon() +   # 比对连接带
  geom_gene() +     # 基因注释
  geom_axis()       # 位置坐标轴
```

![使用全部默认参数的弦图](reference/figures/combined_default.png)

使用全部默认参数的弦图

核心思想就一句话：**数据交给
[`ggchord()`](https://dangjem.github.io/ggchord/reference/ggchord.md)，样式交给各图层**。

## 延伸阅读

| 需求 | 资源 |
|----|----|
| 完整教程：数据准备、导入、校验、清理、示例与灵活参数格式 | [**ggchord 教程**](https://dangjem.github.io/ggchord/articles/ggchord_vignette.html)（[源文件](https://dangjem.github.io/ggchord/vignettes/ggchord_vignette.Rmd)） |
| 中文完整使用指南 | [**ggchord 中文指南**](https://dangjem.github.io/ggchord/articles/ggchord_guide_hans.html)（[源文件](https://dangjem.github.io/ggchord/vignettes/ggchord_guide_hans.Rmd)） |
| 完整的函数与参数参考 | [**函数参考**](https://dangjem.github.io/ggchord/reference/index.html) |
| 版本更新记录 | [**NEWS.md**](https://dangjem.github.io/ggchord/NEWS.md) |

## 贡献与反馈

欢迎通过 [GitHub Issues](https://github.com/DangJem/ggchord/issues) 提交
Bug 报告与功能建议，也欢迎提交 Pull Request。

## 许可证

本项目基于 MIT 许可证发布，详见
[LICENSE](https://dangjem.github.io/ggchord/LICENSE) 文件。
