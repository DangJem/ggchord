先把四条随包序列排成横行。这里用与弦图相同的输入，先观察骨架。

## 绘制四条骨架 {#draw}

按[安装页](installation.html)安装本网站对应版本后，在全新 R 会话运行本课。

```r
library(ggchord)
```

{{code:frontlinearbase}}

{{figure:frontlinearbase}}

## 按共同尺度读位置和长度 {#read}

每一行代表一条序列。比较 `MT108731.1` 和 `MT118296.1`，前者约 64 kb，后者约 32 kb。共线性图用共同的 bp 尺度，因此可以直接比较横向长度。

## 试着去掉序列名 {#try}

去掉 `geom_seq_label()` 后，四条骨架及其相对长度仍在。

## 换用自己的序列 {#adapt}

提供含 `accver` 和完整 `length` 的序列表。接下来[添加比对和 CDS](collinear-build.html)。
