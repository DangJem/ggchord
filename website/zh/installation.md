<div class="lesson-intro">本页完成安装和环境确认，然后选择一个模块逐步学习；不要求先读完全部功能参考。</div>

## 安装并确认环境 {#requirements}
需要R ≥4.1.0、ggplot2 ≥4.0.0。安装和示例均对应v0.16.6。

```r
install.packages(c("remotes", "dplyr", "purrr", "tidyr"))
remotes::install_github("DangJem/ggchord@v0.16.6")
library(ggchord)
packageVersion("ggchord")
```

### 你正在使用哪个版本？ {#source}
基础安装和示例对应 **ggchord 0.16.6**。上面的安装命令固定到该版本的正式标签；发布资产可在 [GitHub Releases](https://github.com/DangJem/ggchord/releases/tag/v0.16.6) 查阅。
如果你使用的是源码 checkout，在仓库根目录运行：

```r
install.packages("pkgload")
pkgload::load_all(".")
packageVersion("ggchord")
```

选择与你手中源码相符的加载方式，不需要同时运行两条安装路径。复现图时记录包版本。
本教程不需要可选的 Bioconductor 或 Biopython 依赖。

### 支持证据 {#support-evidence}
声明下限为 R 4.1.0、ggplot2 4.0.0。发布检查及实际平台版本记录在正式验收中，不代表所有依赖组合均获认证。可选读取器需要其文档列出的依赖。

## 选择一条学习路线 {#first-steps}
1. [模块一『弦图』：先认识三张表](chord-data.html)，再读导入、骨架、图层与参数。
2. [模块二『共线性图』：复用弦图数据](collinear-data.html)，再替换坐标、路由与局部比较。
3. [模块三『质粒图』：从单序列与DNA开始](plasmid-data.html)，再逐层添加，遍历13个参考案例。

建议从弦图开始，符合ggchord的发展顺序。完整学习地图见[全部教程](articles.html)。下面的两序列练习仅用于快速确认环境，不代替三个模块的递进课程。

## 输入表先看哪些列？ {#input}
序列表最少accver和length；连接表最少两个ID及四个端点；基因表最少accver/start/end/strand。[模块一的第一课](chord-data.html)提供真实表格、可选列与BLAST命名说明。

## 可选：快速确认一次绘图 {#first-plot}
运行整个折叠块，再查看图片。原始记录来自随包数据，不需要BLAST或外部文件。

{{prepare:firstdata,firstplot}}
{{figure:first}}

## 如何阅读这个快速结果 {#read}
箭头保留真实链方向，颜色只区分所属序列；连线是已有比对区间，绘图没有生成新的同源结论。该小图选择线性视图便于读位置；完整弦图教学从序列骨架开始。

## 可选：保存快速结果 {#save}
仅当上方代码已运行、p_first已创建时执行：

{{code:firstsave}}

预览与保存使用相同10×5 in尺寸，参见[输出](output.html)。

## 常见运行问题 {#troubleshooting}
找不到对象时，检查是否运行了本页准备；找不到函数时先确认版本。保存位置由getwd()决定。其他格式的可选依赖只在相应专题需要，基础模块不要求安装Bioconductor或Biopython。
