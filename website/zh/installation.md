<a id="requirements"></a>

## 安装 ggchord {#requirements}

已安装 ggchord？[用随包数据试画](#first-plot)。

ggchord 需要 R ≥ 4.1.0 和 ggplot2 ≥ 4.0.0。使用开发预览时，从生成该预览的源码安装。使用已发布版本时，在 [GitHub Releases](https://github.com/DangJem/ggchord/releases) 选择与文档对应的版本，下载并解压源码。

将 R 的工作目录设为包含 `DESCRIPTION` 的 ggchord 源码目录，然后运行：

```{.r}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_local(".", dependencies = NA, upgrade = "never", build = FALSE)
```

已有早期版本的脚本可查看[接口迁移](migration.md)。

<a id="first-plot"></a>

## 用随包数据画一张弦图 {#first-plot}

复制下面整段到 R 中运行。示例使用随包数据，无需另外下载文件。

```r
library(ggchord)

p <- ggchord(data = seq_data_example) +
  geom_seq() +
  geom_seq_label() +
  coord_chord()
p
```

你会看到四条带名称的弧，每条弧对应一条序列，弧长反映完整序列长度。用 `geom_seq()` 绘制骨架，用 `geom_seq_label()` 加上序列名，并用 `coord_chord()` 选择弧形布局。

<a id="first-steps"></a>

## 按你的数据继续 {#first-steps}

- [添加连接和注释](multi-sequence-start.md#links)
- [选择另一种图](modes.md#modules)
- [使用自己的文件](formats.md#choose-format)
- [预览与保存](output.md#size)

<a id="environment"></a>

## 查看当前安装版本 {#environment}

```r
getRversion()
packageVersion("ggplot2")
packageVersion("ggchord")
find.package("ggchord")
```

<a id="support-evidence"></a>

## 示例需要哪些额外依赖？ {#support-evidence}

额外的 R 包或文件解码器会在对应示例中声明。需要读取文件时，先查[支持的格式](formats.md#choose-format)。

<a id="troubleshooting"></a>

## 安装或绘图失败时 {#troubleshooting}

先看错误中点名的依赖是否可用，以及 R 加载的 ggchord 是否来自预期库路径。源码安装报错时保留第一条编译错误。绘图的警告和报错可按[诊断与排错](diagnostics.md)检查。

<span id="input"></span><span id="read"></span><span id="save"></span><span id="source"></span><span id="stable"></span><span id="verified-versions"></span>
