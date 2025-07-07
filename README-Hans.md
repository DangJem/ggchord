🌐 语言切换: 【[现代汉语（Han）](README-Hans.md) | [英文（English）](README.md)】

# ggchord：多序列BLAST比对弦图可视化工具

## 概述
`ggchord` 是一个基于 `ggplot2` 的R函数，用于将多序列的BLAST比对结果可视化为直观的弦图，支持丰富的样式定制，可轻松展示序列间的同源区域与结构关系。`ggchord` 0.1.0版本实现了从简单多序列弦图到**功能更加丰富的**多序列弦图的突破性升级，能够同时展示多条序列间的比对关系：
- 每条序列以圆弧或定制轨迹呈现，按比例映射长度。
- 彩色连接带（ribbon）表示序列间的比对区域，支持按相似度或来源着色。
- 配备可定制的坐标轴，精准标注序列位置与长度。
- 支持整体旋转、序列方向调整等布局优化，适配不同分析场景。

适用于比较基因组学、泛基因组分析、噬菌体 - 宿主序列关系等研究，帮助快速挖掘序列间的同源模式。

## 主要功能
- **多序列支持**：同时展示2条及以上序列的比对关系，无需局限于双序列对比。
- **序列级定制**：
  - 自定义序列顺序、方向（正向/反向）、间隙与半径。
  - 自动或手动指定序列颜色与标签，提升可读性。
- **精细化坐标轴**：
  - 每条序列配备独立坐标轴，含主/次刻度，清晰标注长度位置。
  - 支持调整刻度长度、标签大小与偏移量，平衡美观与信息密度。
- **灵活的连接带样式**：
  - 3种配色方案（单一颜色、按查询序列、按相似度渐变）。
  - 可调节连接带与序列的间隙，支持贝塞尔曲线控制点定制平滑度。
- **布局优化**：整体图形可旋转，适配不同展示需求。
- **调试模式**：辅助排查数据问题，显示有效/无效比对数量。

## 安装
### 依赖环境
- R (≥ 3.6.0)
- ggplot2 (≥ 3.3.0)
- ggnewscale (≥ 0.5.0)
- RColorBrewer

```r
install.packages("ggplot2")
install.packages("ggnewscale")
install.packages("RColorBrewer")
```

### 如何安装ggchord？
从CRAN安装gggenes的稳定版本：

`install.packages("ggchord")`

如果你想要开发版本，请从GitHub安装：

`devtools::install_github("DangJem/ggchord")`

## 使用说明
### 前期数据准备
需准备两类输入数据：

#### 【必须】序列信息数据（`seq_data`）
TSV（制表符分隔值）文件，包含序列基本信息，必须包含以下列：
- `seq_id`：序列唯一标识（如基因名、登录号）
- `length`：序列长度（正数）

示例：
```r
seq_data <- read.delim("seq_track.tsv", sep = "\t", stringsAsFactors = FALSE)
```
其中，`seq_track.tsv` 文件格式如下（示例）：
```txt
seq_id	length
MT108731.1	64323
MT118296.1	32090
OQ646790.1	57367
OR222515.1	83080
```
你可以使用以下命令从FASTA文件自动生成该表格：
```bash
seqkit fx2tab -nil *fna | sed '1i seq_id\tlength' > seq_track.tsv
```

#### 【可选】比对数据（`ribbon_data`）
TSV（制表符分隔值）文件，包含BLAST比对结果（可通过`outfmt6`或`outfmt7`格式转换），必须包含以下列：
- `qaccver`：查询序列ID（需存在于`seq_data$seq_id`）
- `saccver`：目标序列ID（需存在于`seq_data$seq_id`）
- `length`：比对长度
- `pident`：序列相似度（百分比）
- `qstart`/`qend`：查询序列上的比对起止位置
- `sstart`/`send`：目标序列上的比对起止位置

可以用以下脚本对示例序列进行BLAST比对，获得outfmt7格式的结果文件：
```bash
# 使用示例样本的fasta文件进行blast比对的脚本
seqs=("MT108731.1" "MT118296.1" "OQ646790.1" "OR222515.1")
seqsNum=${#seqs[@]}
ext="fna"
for ((i=0; i<seqsNum-1; i++)); do
  for ((j=i+1; j<seqsNum; j++)); do
    echo -e "Running BLASTN: ${seqs[$i]} vs ${seqs[$j]}"
    blastn \
      -outfmt '7 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen sstrand stitle' \
      -query "${seqs[$i]}.${ext}" \
      -subject "${seqs[$j]}.${ext}" \
      -out "${seqs[$i]}__${seqs[$j]}.o7"
  done
done
```

#### 【可选】基因数据（`gene_data`）
TSV（制表符分隔值）文件，该数据框必须包含以下列：
- `seq_id`：序列唯一标识，需与序列信息数据（`seq_data`）中的 `seq_id` 对应，例如基因名、登录号等。
- `start`：基因在序列上的起始位置。
- `end`：基因在序列上的结束位置。
- `strand`：基因所在的链方向，通常用 `+` 或 `-` 表示正向或反向。
- `anno`：基因的注释信息，如基因的功能描述。

示例文件 `gene_track.tsv` 的格式如下：
```txt
seq_id	start	end	strand	anno
MT108731.1	100	200	+	DNA binding protein
MT118296.1	300	400	-	Transcription factor
```

你可以通过运行 `gff2gene_track.R` 脚本，将 gff3 格式的文件转换成符合要求的基因数据表格。以下是脚本的具体内容：
```r
library(tidyverse)

# 获取当前目录下所有 gff3 文件的路径
gff3FilesPath <- list.files(path = ".", pattern = "*.gff3")

# 读取所有 gff3 文件并合并成一个数据框
gff3Table <- map_df(gff3FilesPath,~read_tsv(.x,show_col_types = F,comment = "#",col_names = F) %>% set_names(c("seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes")))

# 筛选出 CDS 类型的记录，并提取注释信息
geneTrackTable <- gff3Table %>% filter(type=="CDS") %>% mutate(anno=str_extract(attributes,"(?<=product=)[^;]+(?=;)")) %>% select(seq_id,start,end,strand,anno)

# 将处理后的数据框保存为 TSV 文件
write_tsv(geneTrackTable,"gene_track.tsv")
```
运行上述脚本后，会在当前目录下生成 `gene_track.tsv` 文件，该文件即可作为基因数据用于后续的分析和可视化。


## 使用示例
### 数据读取
```R
#读取序列长度数据
seq_data <- read.delim("seq_track.tsv", sep = "\t", stringsAsFactors = FALSE)

# 读取并处理BLAST数据
read_blast <- function(file) {
  df <- read.delim(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
  colnames(df) <- c("qaccver","saccver","pident","length","mismatches","gapopen",
                    "qstart","qend","sstart","send","evalue","bitscore",
                    "qcovs","qlen","slen","sstrand","stitle")
  df
}
blast_files <- list.files(path = ".", pattern = "*.o7", full.names = TRUE)
all_blast <- do.call(rbind, lapply(blast_files, read_blast))
ribbon_data <- subset(all_blast, length >= 100)

# 读取基因注释数据，为了使图像更加美观，这里过滤掉了较短的基因注释
gene_data <- read.delim("gene_track.tsv", sep = "\t", stringsAsFactors = FALSE) |> dplyr::slice_max(order_by = end-start, n = 5, by = seq_id)
```

### 只传入必要的seq_data
对于ggchord来说，序列数据是最重要的、不可或缺的。默认情况下，序列将按你所传入的seq_data顺序，以逆时针来排列。当然这些都是可以修改的
```
part1_1 <- ggchord(
  seq_data = seq_data,
)
```R
![plot](/examples/plots/v1.0.0/part1_1.jpg)


例如，在下面的例子中，你可以通过seq_order、seq_orientation及seq_curvature来控制序列的顺序、方向和曲率，通过seq_colors来设置序列的颜色
```R
part1_2 <- ggchord(
  seq_data = seq_data,
  seq_order = c("MT118296.1", "OR222515.1", "MT108731.1", "OQ646790.1"),
  seq_orientation = c(1,-1,1,-1),
  seq_curvature = c(0,2,-2,6),
  seq_colors = c("steelblue", "orange", "pink", "yellow")
)
```
![plot](/examples/plots/v1.0.0/part1_2.jpg)

### 加入序列比对数据
对于基因比对弦图来说，序列比对绝对是我们主要关注的点，因此ribbon_data是除seq_data之外最重要的数据。
默认情况下，ribbon的填充色由blast结果中的百分比一致性决定。
```R
part2_1 <- ggchord(
  seq_data = seq_data,
  ribbon_data = ribbon_data
)
```
![plot](/examples/plots/v1.0.0/part2_1.jpg)

当然，这些也都是可以更改的。例如，你可以将填充色制定为从查询序列出发的颜色，这样用户可以更加容易捕获到不同序列间的比对
```R
part2_2 <- ggchord(
  seq_data = seq_data,
  ribbon_data = ribbon_data,
  ribbon_color_scheme = "query"
)
```
![plot](/examples/plots/v1.0.0/part2_2.jpg)

如果你觉得颜色无所谓，也可以设置成单色
```R
part2_3 <- ggchord(
  seq_data = seq_data,
  ribbon_data = ribbon_data,
  ribbon_color_scheme = "single",
  ribbon_colors = "orange"
)
```
![plot](/examples/plots/v1.0.0/part2_3.jpg)


另外，ribbon会根据序列的方向、弯曲度、间距及半径等参数自动调整以完美匹配（注：坐标轴和基因箭头亦是如此）
> 当前版本仍存在一些问题，再某些参数组合时可能出现图像变形的情况，这将会在未来的版本中得到修复
```R
part2_4 <- ggchord(
  seq_data = seq_data,
  ribbon_data = ribbon_data,
  seq_orientation = c(1,-1,1,-1),
  seq_curvature = c(0,2,-2,6),
  seq_gap = c(.1,.05,.09,.05),
  seq_radius = c(1,5,1,1)
)
```
![plot](/examples/plots/v1.0.0/part2_4.jpg)


## 加入基因注释信息
通过基因注释信息，我们可以更加容易得找到序列的基因区域的一致性区域，以挖掘出不同物种之间的进化关系
```R
part3_1 <- ggchord(
  seq_data = seq_data,
  gene_data = gene_data
)
```
![plot](/examples/plots/v1.0.0/part3_1.jpg)

默认情况下箭头的填充色为链模式，即使用其所在基因序列的链方向进行填充，也可以使用手动指定模式，此时将使用gene_data中的anno注释类别进行颜色填充。
```R
part3_2 <- ggchord(
  seq_data = seq_data,
  gene_data = gene_data,
  gene_color_scheme = "manual"
)
```
![plot](/examples/plots/v1.0.0/part3_2.jpg)

当然，基因注释标签也是可以显示的。不过，如你所间，要想调节出一个完美的效果或许会花费一番功夫
```R
part3_3 <- ggchord(
  seq_data = seq_data,
  gene_data = gene_data,
  gene_label_show = T,
  gene_label_rotation = 45,
  gene_label_radial_offset = .1,
  panel_margin = list(l=.2)
)
```
![plot](/examples/plots/v1.0.0/part3_3.jpg)


## 综合示例
我们一般绘图会同时使用seq_data、ribbon_data和gene_data，这样的图像才会更加赏心悦目。
```R
part4_1 <- ggchord(
  seq_data = seq_data,
  ribbon_data = ribbon_data,
  gene_data = gene_data,
)
```
![plot](/examples/plots/v1.0.0/part4_1.jpg)

当然，ggchord也提供了比较丰富的参数来控制图像的各种细节。
```sh
part4_2 <- ggchord(
  seq_data = seq_data,
  ribbon_data = ribbon_data,
  gene_data = gene_data,
  title = "Multi-sequence Chord Diagram with Gene Annotations",
  seq_gap = .03,
  seq_radius = c(3,2,2,1),
  seq_orientation = c(-1, -1, -1, 1),
  seq_curvature = c(0,1,-1,1.5),
  gene_offset = list(c("+"=.2,"-"=-.2), 
                     c("+"=.2,"-"=-.2), 
                     c("+"=.2,"-"=0), 
                     c("+"=.2,"-"=.1)),
  gene_label_rotation = list(c("+"=45,"-"=-45), 
                             c("+"=.2,"-"=-.2), 
                             c("+"=.2,"-"=0), 
                             c("+"=.2,"-"=.1)),
  gene_label_radial_offset = c(0,0,0,0),
  gene_label_circum_offset = c(1, 0, -2, 0),
  gene_label_circum_limit = c(T,T,T,T),
  gene_width = .08,
  gene_label_show = T,
  gene_color_scheme = "strand",
  ribbon_gap = .1,
  ribbon_color_scheme = "pident",
  ribbon_ctrl_point = c(0,0),
  axis_label_orientation = c(0,45,80,130),
  axis_gap = 0,
  axis_tick_major_number = 5,
  axis_tick_major_length = 0.03,
  axis_tick_minor_number = 5,
  axis_tick_minor_length = 0.01,
  axis_label_size = 2,
  axis_label_offset = 2,
  rotation = 45, 
  show_axis = T,
  panel_margin = list(t=.1),
  debug = TRUE,
)
```
![plot](/examples/plots/v1.0.0/part4_2.jpg)

不过，当前还是软件的早期版本，很多地方还有不完善的地方甚至是BUG。随着日后软件迭代，相信这些问题都将一一解决。


## 参数详情
| 参数类别 | 参数名 | 类型 | 默认值 | 描述 |
| --- | --- | --- | --- | --- |
| **核心数据** | `seq_data` | data.frame/tibble | - | 包含序列信息的数据框，必须包含列：<br> - `seq_id`：序列唯一标识<br> - `length`：序列长度 |
|  | `ribbon_data` | data.frame/tibble | - | 包含BLAST比对结果的数据框，可选。若提供，则必须包含列：<br> - `qaccver`：查询序列ID<br> - `saccver`：目标序列ID<br> - `length`：比对长度<br> - `pident`：序列相似度百分比<br> - `qstart`：查询序列起始位置<br> - `qend`：查询序列结束位置<br> - `sstart`：目标序列起始位置<br> - `send`：目标序列结束位置 |
|  | `gene_data` | data.frame/tibble | - | 包含基因注释信息的数据框，可选。若提供，则必须包含列：<br> - `seq_id`：序列唯一标识<br> - `start`：基因起始位置<br> - `end`：基因结束位置<br> - `strand`：链方向（+或-）<br> - `anno`：基因注释 |
| **基本样式** | `title` | 字符串 | NULL | 图形主标题 |
| **序列布局** | `seq_order` | 字符向量 | NULL | 指定序列绘制顺序；若为NULL则使用 `seq_data` 中的顺序 |
|  | `seq_labels` | 字符向量或命名向量 | NULL | 序列标签；若为NULL则使用 `seq_id` |
|  | `seq_orientation` | 数值向量或单值 | 1 | 每条序列的方向：1（正向）或 -1（反向）；默认正向 |
|  | `seq_gap` | 数值或向量 | 0.03 | 长度与序列数一致，定义每条序列头部到下一条序列尾部的弧度比例 [0,0.5) |
|  | `seq_radius` | 数值或向量 | 1.0 | 序列圆弧半径，支持单值或与序列数相同的向量 |
|  | `seq_curvature` | 数值或向量 | 1.0 | 序列圆弧的弯曲程度：1为标准圆弧，0为直线，>1更弯曲 |
|  | `seq_colors` | 颜色向量或命名向量 | NULL | 定义各序列圆弧颜色；若为NULL则基于 RColorBrewer Set1 自动生成 |
| **基因样式** | `gene_offset` | 数值、向量或列表 | 0.03 | 基因箭头与序列圆弧之间的径向偏移距离。支持：<br> - 单值：所有序列的所有链使用相同偏移<br> - 向量：长度与序列数一致，每个序列的所有链使用相同偏移<br> - 列表：命名列表，每个元素对应一个序列，元素可为单值（该序列所有链）或包含"+"和"-"的命名向量（区分链） |
|  | `gene_width` | 数值或向量 | 0.1 | 基因箭头宽度 |
|  | `gene_label_show` | 逻辑值 | FALSE | 是否显示基因标签 |
|  | `gene_label_rotation` | 数值、向量或列表 | 0 | 基因标签的旋转角度（度），支持与 `gene_offset` 相同的参数格式 |
|  | `gene_label_size` | 数值 | 2.5 | 基因注释文字大小 |
|  | `gene_label_radial_offset` | 数值、向量或列表 | 0 | 基因标签相对于箭头的径向偏移（正值向外，负值向内），支持与 `gene_offset` 相同的参数格式 |
|  | `gene_label_circum_offset` | 数值、向量或列表 | 0 | 基因标签沿序列周向的偏移比例（相对于基因长度），支持与 `gene_offset` 相同的参数格式 |
|  | `gene_label_circum_limit` | 逻辑值、向量或列表 | TRUE | 是否限制周向偏移不超过基因长度的一半，支持与 `gene_offset` 相同的参数格式 |
|  | `gene_color_scheme` | 字符 | "strand" | 指定基因颜色方案，可选"strand"（按链方向）或"manual"（手动指定） |
|  | `gene_colors` | 颜色向量 | - | 用于指定基因箭头的填充色，具体行为取决于 `gene_color_scheme`：<br> - "strand"模式：支持命名向量（仅"+"/"-"）、非命名向量（先"+"后"-"）或单值（正负链同色），缺省则"+"为红色、"-"为蓝色<br> - "manual"模式：支持命名向量（对应 `anno`）、非命名向量（截断多余，补齐不足），缺省则使用自动生成颜色 |
|  | `gene_order` | 字符向量 | NULL | 指定基因在图例中的显示顺序；若为NULL则使用基因在数据中出现的顺序 |
| **连接带样式** | `ribbon_color_scheme` | 字符 | "pident" | 连接带配色方案，可选 "single"、"query" 或 "pident" |
|  | `ribbon_colors` | - | - | 连接带颜色参数：<br> - single：单一颜色（单值或向量取第一个）<br> - query：按查询序列映射颜色（命名或非命名向量或单值）<br> - pident：渐变色阶向量，用于按相似度百分比生成渐变，默认使用蓝到黄的渐变色 |
|  | `ribbon_alpha` | 数值 | 0.35 | 连接带透明度 [0,1] |
|  | `ribbon_ctrl_point` | 向量或列表 | - | 贝塞尔控制点，用于调整连接带形状：<br> - 向量：长度为2（单控制点）或4（c1x,c1y,c2x,c2y，双控制点）<br> - 列表：每个元素为包含1 - 2个控制点的子列表，默认自动计算 |
|  | `ribbon_gap` | 数值或向量 | 0.15 | 序列圆弧与连接带之间的径向距离 |
| **坐标轴设置** | `axis_gap` | 数值或向量 | 0.04 | 坐标轴与序列圆弧之间径向距离，支持负值 |
|  | `axis_tick_major_number` | 整数或向量 | 5 | 每条序列主刻度数 |
|  | `axis_tick_major_length` | 数值或向量 | 0.02 | 主刻度线长度比例 |
|  | `axis_tick_minor_number` | 整数或向量 | 4 | 每两个主刻度之间的次刻度数 |
|  | `axis_tick_minor_length` | 数值或向量 | 0.01 | 次刻度线长度比例 |
|  | `axis_label_size` | 数值或向量 | 3 | 坐标轴刻度文字大小 |
|  | `axis_label_offset` | 数值或向量 | 0 | 坐标轴标签相对于刻度线的偏移比例 |
|  | `axis_label_orientation` | 字符、数值或向量 | "horizontal" | 坐标轴标签方向：<br> - "horizontal"：水平方向<br> - 数值：旋转角度（度）<br> - 向量：长度与序列数一致或命名向量 |
| **布局设置** | `rotation` | 数值 | 45 | 整个图形的旋转角度（度） |
|  | `panel_margin` | 列表 | list(t=0,r=0,b=0,l=0) | 图形边缘留白（t=上, r=右, b=下, l=左） |
| **图例与调试** | `show_legend` | 逻辑值 | TRUE | 是否显示图例 |
|  | `show_axis` | 逻辑值 | TRUE | 是否显示坐标轴及刻度 |
|  | `debug` | 逻辑值 | FALSE | 是否输出调试信息 | 

## 图形解读
- **序列圆弧**：每条彩色圆弧代表一条序列，长度按比例映射，箭头指示方向（正向/反向）。
- **连接带**：连接不同序列的彩色区域，代表比对区间：
  - 按相似度着色时，颜色渐变反映序列一致性（如从蓝色到红色表示相似度从低到高）。
  - 按查询序列着色时，同一颜色的连接带来自同一条查询序列。
- **坐标轴**：每条序列外侧的刻度线与数字，标注序列位置（单位与序列长度一致），方便定位比对区域。

## 版本更新记录
### v0.2.0（最新）
- **高级弧线和直线模式优化**：
  - 通过增强的曲线拟合算法，在双序列可视化中，“弧线模式”和“直线模式”在不同比对区域之间实现了更平滑的过渡。例如，在可视化两条密切相关且具有多个短比对的序列时，弧线或直线现在能够更优雅地连接这些区域，减少视觉混乱。
  - 确保了比对数据的更准确表示。在以前的版本中，尤其是对于长距离比对，可视化可能存在轻微的失真。现在，弧线和直线的长度和位置更精确地匹配实际比对坐标，提供了更真实的序列关系视图。
- **精确的曲率和间隙控制**：
  - 用户现在可以更精确地控制弦图中弧线的曲率。引入了一个新参数，允许逐步调整弧线曲率，使用户能够突出不同类型的比对模式。例如，可以使用更明显的曲率来强调高度保守的区域，而较平坦的曲线则用于不太重要的比对。
  - 序列和连接带之间的间隙可以进行更精细的调整。这对于可视化复杂的比对场景非常有用，其中不同序列具有不同程度的相似性。用户现在可以为不同的序列对设置不同的间隙值，确保可视化既清晰又有信息量。
- **增强的颜色定制**：
  - 0.2.0版本为序列和连接带提供了更广泛的调色板。除了现有的默认配色方案外，用户现在可以从各种预定义的调色板中进行选择，这些调色板针对不同类型的数据可视化进行了优化。例如，有专门设计用于突出高对比度区域的调色板，以及用于创建更微妙和美观可视化的调色板。
  - 在使用`pident`配色方案为连接带着色时，用户现在可以自定义渐变颜色。这允许进行更个性化的可视化，更好地表示序列相似性的分布。例如，用户可以选择使用类似热图的渐变来清晰显示从低到高的相似性值范围。
### v0.1.0
- 支持通过`seq_data`、`ribbon_data`与`gene_data`分离管理序列、比对于基因数据。
- 新增序列方向控制、自定义顺序、间隙与半径调整。
- 实现可定制坐标轴（主/次刻度、标签位置）。
- 连接带支持3种配色方案（单一色、按查询序列、按相似度渐变）。
- 支持整体图形旋转与调试模式。
### v0.0.2
- 支持多序列展示（从双序列升级）
- 新增双序列的“弧线模式”与“直线模式”切换。
- 支持弧线曲率调整与间隙控制。
### v0.0.1
- 初始版本，支持双序列BLAST比对弦图可视化。

## 贡献与反馈
欢迎提交issue报告bug或建议新功能，也可通过Pull Request参与开发。使用中遇到问题，可参考示例代码或开启调试模式（`debug = TRUE`）排查数据问题。

## 许可证
本项目采用MIT许可证，详情见[LICENSE](LICENSE)文件。
