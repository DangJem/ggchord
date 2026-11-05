## 下载与安装 {#requirements}

需要 R ≥ 4.1.0 和 ggplot2 ≥ 4.0.0。在 R 中运行下列命令，CRAN 会提供 ggchord 和所需依赖。

```r
install.packages("ggchord")
library(ggchord)
packageVersion("ggchord")
```

最后一行显示已安装的版本。本教程使用 0.16.6 的公开接口。

### 源码下载 {#source}

可从 [CRAN](https://cran.r-project.org/web/packages/ggchord/index.html) 或 [GitHub v0.16.6](https://github.com/DangJem/ggchord/releases/tag/v0.16.6) 下载源码包。

### 可选依赖 {#support-evidence}

开始使用和基础弦图章节只需要 ggchord。完整案例会列出额外依赖。

## 确认安装 {#troubleshooting}

在新 R 会话中用 `library(ggchord)` 加载。如果安装或加载报错，先查看报错中提到的 R 或依赖包版本。

## 下一步 {#first-steps}

- [基本数据表](chord-data.html)展示序列、连接和特征的真实记录，并逐列解释。
- [绘图基本概念](modes.html)介绍骨架、图层与坐标系。
- <span id="first-plot"></span><span id="input"></span><span id="read"></span>[绘制序列骨架](chord-start.html)进入弦图章节。
- <span id="save"></span>[保存图片](output.html)介绍尺寸与输出格式。
