ggchord 根据序列坐标展示特征的位置，以及不同序列之间的比对关系。假设你有两条噬菌体基因组及它们之间的一条局部比对。可以先按全长画出两条基因组，再标出相互匹配的区间，最后查看这些区间内有哪些编码序列。也可以把同一对基因组平行排开，从左到右比较它们的位置。

<a id="requirements"></a>

## 安装 ggchord {#requirements}

要运行本站示例，请[下载与手册配套的包源码](../../)并解压。在含有 `DESCRIPTION` 的 ggchord 目录中打开 R。使用 R 4.1.0 或更新版本和 ggplot2 4.0.0 或更新版本。示例使用 dplyr 筛选和整理记录。

```{.r}
install.packages(c("ggplot2", "remotes", "dplyr"))
remotes::install_local(".", dependencies = NA, upgrade = "never", build = FALSE)
```

如需安装 CRAN 的稳定发布版，可使用以下命令。运行本手册的示例请使用上面提供的配套源码。

```{.r}
install.packages("ggchord")
```

<a id="layers"></a>

## 分三步读懂序列比较图 {#layers}

第一幅图画出两条完整基因组。长弧是 MT108731.1（64,323 bp），短弧是 MT118296.1（32,090 bp）。第二幅图增加一条连接带，带的两端标出第一条局部比对所报告的两个区间。第三幅图中的绿色箭头标出各条编码序列的注释位置。沿着骨架观察箭头的指向，可以读出编码序列所在的链；顺着连接带看向另一端，可以找到对应区间内的特征。

{{preview:first-backbones}}

{{preview:first-chord}}

{{preview:first-features}}

<a id="layouts"></a>

## 在弧线和平行直线上对照同一条比对 {#layouts}

在弦图上，可以从连接带的任一端沿圆周查看对应的序列区间。在共线性图上，可以沿直线骨架比较该区间与基因组两端的距离。两幅图所用的基因组、比对记录和端点坐标相同。

{{preview:views-chord}}

{{preview:views-collinear}}

[用随包数据绘制第一张弦图](chord.md#backbones)。准备自己的表格时，可先读[数据指南](data.md)；要挑选其他图形，可查阅[绘图指南](modes.md)。

把注释类别映射到填色，或把比对得分映射到颜色，标度就会提供解读这些值的图例。[定量轨道](tracks.md)可以沿骨架展示 GC 含量；[Fusion 组合图](fusion.md)可以把总览中的一个区间连到单独的细节面板。

<a id="first-plot"></a>

<a id="example-first-plot"></a>

要直接运行一条比对的示例，可前往[连接一对比对区间](chord.md#draw)。
