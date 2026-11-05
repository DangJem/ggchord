## 下载与安装 {#requirements}

需要 R ≥ 4.1.0 和 ggplot2 ≥ 4.0.0。在 R 中运行下列命令，CRAN 会提供 ggchord 和所需依赖。

```r
install.packages("ggchord")
library(ggchord)
packageVersion("ggchord")
```

最后一行显示已安装版本。**本网站面向开发版 0.16.7.9000**，CRAN 正式版为 0.16.6（2026-11-05 核查）。两者不能视作完全相同：开发版新增蛋白教学对象、跨行连接专用线宽，并修改圆图零旋转基准和部分默认外观。

开发版尚未推送或发布，请从本网站对应的本地源码 checkout 安装；不要假定远端 main 已包含这些修改。在仓库根目录执行：

```r
install.packages(".", repos = NULL, type = "source")
library(ggchord)
stopifnot(as.character(packageVersion("ggchord")) == "0.16.7.9000")
```

先安装依赖，再在新的 R 会话运行以上代码。版本号不能识别同版本不同提交；构建清单 `build-manifest.json` 和审计记录保留源码 SHA。只有 CRAN 版的读者应使用 [0.16.6 对应离线教程与参考](https://stat.ethz.ch/CRAN/web/packages/ggchord/index.html)，不要运行开发版蛋白数据案例或 `nonadjacent_linewidth` 参数。

### 已实测的版本差异 {#verified-versions}

三种首课的基本绘图代码分别在 CRAN 0.16.6 和当前开发 checkout 的全新 R 进程渲染，共六次。这验证基础调用，不保证外观相同，也不代表后续开发版课程都适用于正式版。以下能力也已针对两个命名空间实测。

| 能力 | 0.16.6 | 0.16.7.9000 |
| --- | --- | --- |
| `geom_link_ribbon(nonadjacent_linewidth=...)` | 不提供 | 提供 |
| 显式 `read_fasta_lengths(sequence_unit=...)` | 不提供 | 提供 |
| 随包原生蛋白数据 | 不提供 | 提供 |

正式版源码包 SHA256 与逐进程结果保存在[版本审计记录](../assets/source/tools/teaching-audit/version-check.json)。

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
