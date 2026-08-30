🌐 语言：[简体中文](README-Hans.md) | [English](README.md)

# ggchord

`ggchord` 是一个基于 `ggplot2` 的 R 包，用于绘制分层的多序列弦图。它可以从
规整的数据表中展示序列弧线、比对连接带、基因组特征、标签和位置坐标轴。

## 安装

```r
install.packages("ggchord")
```

需要开发版本时可从 GitHub 安装：

```r
devtools::install_github("DangJem/ggchord")
```

## 快速开始

```r
library(ggchord)

data(seq_data_example)
data(ribbon_data_example)
data(gene_data_example)

ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
  geom_seq() +
  geom_ribbon() +
  geom_gene() +
  geom_axis()
```

输入数据需要包含以下列：

| 数据 | 必需列 |
| --- | --- |
| 序列 | `seq_id`、`length` |
| 序列比对 | `qaccver`、`saccver`、`length`、`pident`、`qstart`、`qend`、`sstart`、`send` |
| 基因组特征 | `seq_id`、`start`、`end`、`strand`、`anno` |

## 文档

- [中文教程](https://dangjem.github.io/ggchord/articles/ggchord_guide_hans.html)
- [英文教程](https://dangjem.github.io/ggchord/articles/ggchord_vignette.html)
- [函数参考](https://dangjem.github.io/ggchord/reference/index.html)
- [版本记录](NEWS.md)

欢迎通过 [GitHub Issues](https://github.com/DangJem/ggchord/issues) 报告问题或提出建议。

## 许可证

MIT，详见 [LICENSE](LICENSE)。
