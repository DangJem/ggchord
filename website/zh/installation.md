<div class="lesson-intro">
用两条随包序列，完成一次从数据到图片的绘图流程。你不需要先运行 BLAST，也不需要准备自己的文件。
</div>

<div class="lesson-goal">
本页完成后，你会得到一张包含序列骨架、基因箭头和比对区间的共线性图，并把它保存为指定尺寸的 PNG。后面的真实数据教程会沿用这套图层思路。
</div>

<div class="lesson-path">安装与加载 → 检查三张表 → 构建图层 → 阅读结果 → 保存图片</div>

## 1. 安装与加载 {#requirements}
需要 R ≥4.1.0 和 ggplot2 ≥4.0.0。先打开 R 或 RStudio，在控制台运行安装命令。
已经安装并能加载 ggchord 时，可以跳过安装。

```r
install.packages("remotes")
remotes::install_github("DangJem/ggchord@v0.16.5")
library(ggchord)
packageVersion("ggchord")
```

### 你正在使用哪个版本？ {#source}
本教程对应 **ggchord 0.16.5**。上面的安装命令固定到该版本的正式标签；发布资产可在 [GitHub Releases](https://github.com/DangJem/ggchord/releases/tag/v0.16.5) 查阅。
如果你使用的是源码 checkout，在仓库根目录运行：

```r
install.packages("pkgload")
pkgload::load_all(".")
packageVersion("ggchord")
```

选择与你手中源码相符的加载方式，不需要同时运行两条安装路径。复现图时记录包版本。
本教程不需要可选的 Bioconductor 或 Biopython 依赖。

## 2. 看懂绘图需要的三张表 {#input}
我们选择 `MT108731.1` 和 `MT118296.1` 两条序列。完整随包对象保留原始字段，
下面的代码只选本例需要的列。先复制运行整个代码块：

{{code:firstdata}}

| 对象 | 必要列 | 每行表示什么 |
| --- | --- | --- |
| `first_sequences` | `accver`, `length` | 一条序列及其完整长度 |
| `first_genes` | `accver`, `start`, `end`, `strand` | 一条基因注释 |
| `first_links` | 两个 ID 和 query/subject 四个端点 | 两条序列之间的一段已有比对 |

`accver` 是连接三张表的标识，必须一致；不需要一定是 GenBank 编号。坐标使用一基闭区间 bp。
随包基因例子额外保存了 `source_strand`，这里恢复真实链方向，没有为图形交替编造正负链。

运行下面的检查。先看列名，再看每张表的前几行：

```r
names(first_sequences)
names(first_genes)
names(first_links)
head(first_genes)
```

<div class="lesson-check">
检查点：序列表 2 列、基因表 4 列、连接表 6 列。没有 pident 或 qlen 也能完成这张基础图；
只有按相似度配色、计算覆盖率等步骤才需要相应字段。
</div>

## 3. 构建第一张图 {#first-plot}
先用线性视图比较位置，图中的距离更容易直接阅读。
`ggchord()` 接收数据，三个几何层分别绘制连接、骨架和基因；坐标系选择线性排列。
使用参数名可避免误把基因表传入第二个位置参数 `ribbon_data`。

{{code:firstplot}}
{{figure:first}}

两种颜色区分基因所属序列，**不表示它们已经被分入同源基因组**。
本例将连接保留为中性灰色，使基因及位置关系更突出。

## 4. 从图中读出什么？ {#read}
按下面顺序观察：

1. 左侧名称对应两条序列，横向骨架和刻度给出完整长度与 bp 位置。
2. 箭头表示随包筛选出的基因注释及其真实链方向；没有画出的基因不等于不存在。
3. 两行之间的带状连接表示已有比对区间，不是绘图库新计算的比对，也不是逐碱基比对图。

这张图适合学习位置与连接。如果目标是查看多序列之间的总体连接关系，可以改用弦图；
如果目标是一条质粒的注释，使用单序列圆图。先明确要回答的问题，再选择图形，见[三类图怎么选](modes.html)。

## 5. 按明确尺寸保存 {#save}
`p_first` 是普通的 ggplot 对象。下面将它保存到当前工作目录：

{{code:firstsave}}

```r
file.exists("first-comparison.png")
```

返回 `TRUE` 表示文件已经生成。RStudio 的 Plots 面板会随窗口大小改变布局；
需要预览最终排版时，用与导出相同的尺寸：

```r
view_ggchord(p_first, width = 10, height = 5, units = "in")
```

## 6. 解决第一次运行时常见的问题 {#troubleshooting}
| 现象 | 先检查什么 |
| --- | --- |
| 找不到 `first_sequences` 或 `p_first` | 是否从本页第 2 步开始依次运行，而非只复制最后一行 |
| 提示没有某个函数 | 是否加载了与本教程对应的版本；先查看 `packageVersion()` |
| 没有连接 | 两端 ID 是否都出现在序列表，筛选后连接表是否为空 |
| 图很小或文字布局不同 | 先按上面的 10 × 5 in 尺寸预览，见[输出与尺寸](output.html) |
| 保存后找不到图片 | 用 `getwd()` 确认当前工作目录 |

## 接下来完成哪项任务？ {#first-steps}
现在你已经完成一次数据准备、图层组合和保存。按你的下一步需要继续：

- **换成自己的数据**：[读入和转换数据](formats.html)，BLAST 用户继续看[默认表格与蛋白坐标](alignments.html)。
- **比较完整的真实记录**：[FASTA、GFF3 和比对表工作流](real-data.html)。
- **绘制质粒图**：[pUC19 的基因、引物和酶切注释](plasmids.html)。
- **控制外观**：[颜色与图例](styling.html)、[标签安排](labels.html)。
- **查某个参数**：[完整函数参考](reference.html)。
