<a id="first-plot"></a>

## 绘制第一张弦图 {#first-plot}

序列表给出四条噬菌体基因组序列的长度，比对表的前四行给出由连接带相连的区间对。安装 ggchord 后，在 R 中运行这五行代码。

{{example:first-plot}}

<a id="requirements"></a>

## 从源码安装 {#requirements}

网站示例使用本手册配套的源码版本，需要 R 4.1.0 或更高版本以及 ggplot2 4.0.0 或更高版本。在含有 `DESCRIPTION` 的 ggchord 源码目录中打开 R，按需安装 remotes 包，再运行本地安装命令。

```{.r}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_local(".", dependencies = NA, upgrade = "never", build = FALSE)
```

安装后重启 R，用 `packageVersion("ggchord")` 核对版本。要换成自己的序列长度和比对区间对，见[弦图指南](chord.md)。
