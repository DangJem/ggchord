## 从 GitHub 安装 {#requirements}

从项目仓库安装 ggchord 及必需依赖。源码安装需要可用的 R 开发环境。

```r
install.packages("remotes")
remotes::install_github("DangJem/ggchord", dependencies = NA)
```

## 加载与绘制 {#first-plot}

```r
library(ggchord)

p <- ggchord(seq_data_example) +
  geom_seq() +
  geom_seq_label() +
  coord_chord()
p
```

## 可选依赖与输入文件 {#support-evidence}

使用 dplyr、Biostrings 或格式解码器的示例会在代码旁声明依赖，按对应工作流安装即可。随包数据可离线使用；文件示例使用[输入归档](formats.html)。

## 安装失败时 {#troubleshooting}

先检查第一条依赖或编译错误，安装后重新启动 R，再用 `help(package = "ggchord")` 查看可用函数。绘图问题见[诊断与排错](diagnostics.html)。

## 相关说明 {#first-steps}

[绘图基本概念](modes.html) · [输入数据表](data.html) · [保存图形](output.html)


<span id="input"></span><span id="read"></span><span id="save"></span><span id="source"></span><span id="stable"></span><span id="verified-versions"></span>
