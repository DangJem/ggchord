## 下载与安装 {#requirements}

需要 R ≥ 4.1.0 和 ggplot2 ≥ 4.0.0。在 R 中运行下列命令，CRAN 会提供 ggchord 和所需依赖。

```r
install.packages("ggchord")
library(ggchord)
packageVersion("ggchord")
```

最后一行显示已安装版本。本地文档预览 **0.25.0.9000**；仓库记录的最近正式版为
**0.16.6**。新的基础架构、Fusion、BEDPE与字符对齐接口需要相应开发 checkout。更新已有脚本前请阅读
[迁移说明](../assets/MIGRATION-Hans.md)。

开发时先安装依赖，再在对应 checkout 根目录使用新的 R 会话运行：

```r
install.packages(".", repos = NULL, type = "source")
library(ggchord)
stopifnot(as.character(packageVersion("ggchord")) == "0.25.0.9000")
```

版本号不能区分不同提交，下载的构建清单记录源码版本。本地预览不代表包已发布或
网站已部署。正式版用户可查阅 [0.16.6 文档](https://stat.ethz.ch/CRAN/web/packages/ggchord/index.html)。

<span id="verified-versions"></span>

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
