[English](README.md) · [现代汉语](README-Hans.md)

<img src="man/figures/logo.png" alt="ggchord" width="150" align="right">

# ggchord：弦图、圆图与共线性图

基于 ggplot2 的 R 包。用独立图层组合序列骨架、区间连接、基因、通用注释、标签和定量轨道。

这里展示**0.16.5 开发版本（基于 0.16.4 候选 API）**，当前正式发布版仍为 0.15.0。源码中的新网站尚未部署。

## 安装

R ≥ 4.1.0 · ggplot2 ≥ 4.0.0

```r
install.packages("remotes")
remotes::install_github("DangJem/ggchord")
```

如已检出本版源码，在仓库根目录使用 `pkgload::load_all(".")`。不声称存在 CRAN 发布。

![真实序列与全部注释的共线性案例](man/figures/reallinear.png)

图中使用仓库原始 FASTA、GFF3 与 BLAST；完整读取、筛选代码与解释见[真实数据专题](vignettes/ggchord_guide_hans.Rmd#real-data)。


## 三种模式

下列为合成教学输入，不代表生物学观测。

```r
library(ggchord)
seqs <- data.frame(accver = c("A", "B", "C"), length = c(1000, 800, 1200))
links <- data.frame(qaccver = c("A", "B"), saccver = c("B", "C"),
  qstart = c(100, 200), qend = c(350, 500),
  sstart = c(150, 400), send = c(400, 700), pident = c(96, 89))
genes <- data.frame(accver = c("A", "B", "C"),
  start = c(150, 200, 400), end = c(300, 350, 650),
  strand = c("+", "+", "-"), anno = c("target", "target", "target"))
```

```r
p_chord <- ggchord(seqs, links, genes) +
  geom_seq() + geom_link_ribbon() + geom_gene() + geom_seq_label() +
  coord_chord()
```

```r
features <- data.frame(accver = "A", start = c(80, 380, 700),
  end = c(280, 610, 850), strand = c("+", "-", "."),
  type = c("CDS", "CDS", "region"), label = c("alpha", "beta", "origin"))
p_circular <- ggchord(seqs[1, ]) + geom_seq() +
  geom_feature(data = features) +
  geom_feature_label_repel(data = features) + coord_circular() +
  theme_ggchord(axis = element_blank())
```

```r
p_collinear <- ggchord(seqs, links, genes) +
  geom_link_ribbon(alpha = .3) + geom_seq() +
  geom_gene(position = position_strand(.1)) + geom_seq_label() +
  coord_collinear(anchor = genes[c("accver", "start", "end")],
    align = "center", orientation = "auto")
```

| Chord / 弦图 | Circular / 圆图 | Collinear / 共线性图 |
| --- | --- | --- |
| ![Chord](man/figures/chord.png) | ![Circular](man/figures/circular.png) | ![Collinear](man/figures/collinear.png) |

`aes()` 映射数据，`scale_*()` 定义映射和图例，geom 定义几何与固定样式，Position 负责局部放置，coord 负责全局坐标，theme 负责非数据外观。

## 预览与保存

RStudio Plots 会按面板尺寸重新排列标签。按相同显式宽高预览并保存：

```r
view_ggchord(p_collinear, width = 8, height = 5, units = "in")
ggsave("comparison.png", p_collinear, width = 8, height = 5, units = "in", dpi = 144)
```

## 继续阅读

- [完整用户教程](vignettes/ggchord_guide_hans.Rmd)：真实数据、质粒、输入格式、连接、标签、轨道、投影与同源、样式、输出和 FAQ。
- [迁移指南](MIGRATION-Hans.md)与 [NEWS](NEWS.md)。
- [本地网站构建与预览](website/README.md)，含函数分类、术语、图集与部署方案。

[Issues](https://github.com/DangJem/ggchord/issues) · [MIT license](LICENSE)
