## 环境要求 {#requirements}
使用 R ≥ 4.1.0 和 ggplot2 ≥ 4.0.0。ggchord 本身无需编译代码。
可选的区间适配器使用 Bioconductor 包，表格导入使用 readxl；只安装实际输入格式
需要的可选依赖。

## 选择来源 {#source}
本指南使用 GitHub 仓库安装。开发 checkout 可能包含正式版之后的变更，复现时应记录
提交号。本指南不声称包已在 CRAN 发布。

```r
install.packages("remotes")
remotes::install_github("DangJem/ggchord")
library(ggchord)
packageVersion("ggchord")
```

对已检出的 0.16.5 源码，在仓库根目录运行：

```r
install.packages("pkgload")
pkgload::load_all(".")
```

[文件导入示例](data.html#import)中的仓库路径不属于安装后的随包数据。
入门例子直接在 R 中构造小型数据框，安装后即可使用。

## 开始使用 {#first-steps}
先读[数据契约](data.html)，复制[共享准备代码](index.html#quick-start)，再选择
[坐标系](modes.html)。保存时显式指定宽高；调整 RStudio Plots 面板会使标签重新排列。
