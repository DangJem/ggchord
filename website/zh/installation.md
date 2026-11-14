## 从 GitHub 安装 {#requirements}

从项目仓库安装 ggchord 及必需依赖。源码安装需要可用的 R 开发环境。

```r
install.packages("remotes")
remotes::install_github("DangJem/ggchord", dependencies = NA)
```

## 环境与安装检查 {#environment}

需要 R ≥ 4.1.0、ggplot2 ≥ 4.0.0。安装后可查看当前会话实际加载的环境；排查问题时一并提供这些输出。

```r
getRversion()
packageVersion("ggplot2")
packageVersion("ggchord")
find.package("ggchord")
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

成功时应看到四条长度不同、标有序列 ID 的弧形骨架，且没有构建错误。若函数缺失，先检查库路径和会话加载的包，再核对对应[函数参考](reference/index.html)。

## 可选依赖与输入文件 {#support-evidence}

使用 dplyr、Biostrings 或格式解码器的示例会在代码旁声明依赖，按对应工作流安装即可。随包数据可离线使用；文件示例使用[输入归档](formats.html)。

## 安装失败时 {#troubleshooting}

| 症状 | 先检查 | 可操作的处理 |
| --- | --- | --- |
| 无法下载仓库或索引 | 第一条网络错误、代理与 GitHub 可达性 | 恢复网络后重试，不先修改示例数据 |
| 依赖缺失或版本不足 | 错误中点名的依赖及 R 下限 | 安装或更新该依赖，再重试 ggchord |
| 源码编译失败 | 日志中第一条编译错误、平台构建工具 | 补齐对应平台 R 开发工具，保留完整日志 |
| 安装成功但函数缺失 | `find.package()` 与新会话加载路径 | 重启 R，确认加载正确库后重试最小图 |

绘图和输入问题见[诊断与排错](diagnostics.html)。

## 相关说明 {#first-steps}

[绘图基本概念](modes.html) · [输入数据表](data.html) · [保存图形](output.html)


<span id="input"></span><span id="read"></span><span id="save"></span><span id="source"></span><span id="stable"></span><span id="verified-versions"></span>
