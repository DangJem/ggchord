## 选择版本与依赖 {#requirements}

本教程使用 **0.27.0.9000 开发接口**。先选择匹配的安装来源；正式版 0.16.6 与开发版不可混用教学范围。

R ≥ 4.1.0; ggplot2 ≥ 4.0.0.

```r
install.packages("remotes")
```

### 正式版 {#stable}

最近正式版为 0.16.6。其源码与对应版本帮助见 [GitHub v0.16.6](https://github.com/DangJem/ggchord/releases/tag/v0.16.6)。安装固定标签后，请使用该版本附带帮助；本教程含较新的接口。

```r
remotes::install_github("DangJem/ggchord@v0.16.6", dependencies = NA)
```

### 开发快照 {#source}

开发运行时固定为公开审计提交 `0f5dd96da57876da74cd1234c7de9cf5579c80b8`，不要从最大命中位置推断版本。当前本地文档重构没有改变运行时行为；文档提交、构建时间与输入哈希另见[构建清单](../build-manifest.json)。安装需联网下载源码和依赖。

```r
runtime_commit <- "0f5dd96da57876da74cd1234c7de9cf5579c80b8"
remotes::install_github(
  paste0("DangJem/ggchord@", runtime_commit),
  dependencies = NA, upgrade = "never"
)
```

<span id="verified-versions"></span>

## 新会话确认与第一图 {#troubleshooting}

```r
library(ggchord)
stopifnot(as.character(packageVersion("ggchord")) == "0.27.0.9000")
stopifnot(packageDescription("ggchord")$RemoteSha ==
  "0f5dd96da57876da74cd1234c7de9cf5579c80b8")
p <- ggchord(seq_data_example) +
  geom_seq() +
  geom_seq_label() +
  coord_chord()
p
```

### 可选依赖 {#support-evidence}

基础首图只需 ggchord；数据整理课程显式使用 dplyr，各导入课单列可选解码器。

## 继续学习 {#first-steps}

<span id="first-plot"></span><span id="input"></span><span id="read"></span>

- [基本数据表](chord-data.html)
- [第一张多序列图](chord-start.html)
- <span id="save"></span>[按目标尺寸保存](output.html)
