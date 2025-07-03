# ggchord: 多序列BLAST比对弦图可视化工具  

基于ggplot2的R函数，用于将多序列的BLAST比对结果可视化为直观的弦图，支持丰富的样式定制，轻松展示序列间的同源区域与结构关系。  


## 概述  
`ggchord` 0.1.0版本实现了从双序列到**多序列弦图**的突破性升级，能够同时展示多条序列间的比对关系：  
- 每条序列以圆弧或定制轨迹呈现，按比例映射长度  
- 彩色连接带（ribbon）表示序列间的比对区域，支持按相似度或来源着色  
- 配备可定制的坐标轴，精准标注序列位置与长度  
- 支持整体旋转、序列方向调整等布局优化，适配不同分析场景  

适用于比较基因组学、泛基因组分析、噬菌体-宿主序列关系等研究，帮助快速挖掘序列间的同源模式。  


## 主要功能  
- **多序列支持**：同时展示2条及以上序列的比对关系，无需局限于双序列对比  
- **序列级定制**：  
  - 自定义序列顺序、方向（正向/反向）、间隙与半径  
  - 自动或手动指定序列颜色与标签，提升可读性  
- **精细化坐标轴**：  
  - 每条序列配备独立坐标轴，含主/次刻度，清晰标注长度位置  
  - 支持调整刻度长度、标签大小与偏移量，平衡美观与信息密度  
- **灵活的连接带样式**：  
  - 3种配色方案（单一颜色、按查询序列、按相似度渐变）  
  - 可调节连接带与序列的间隙，支持贝塞尔曲线控制点定制平滑度  
- **布局优化**：整体图形可旋转，适配不同展示需求  
- **调试模式**：辅助排查数据问题，显示有效/无效比对数量  


## 安装方法  
### 依赖包  
- R (≥ 3.6.0)  
- ggplot2 (≥ 3.3.0)  
- RColorBrewer (≥ 1.1-3)  
- grDevices (内置包)  

### 安装步骤  
```r
# 安装依赖
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
```


## 使用说明  

### 数据准备  
需准备两类输入数据：  

1. **序列信息数据（seq_data）**  
   数据框，包含序列基本信息，必须包含以下列：  
   - `seq_id`：序列唯一标识（如基因名、登录号）  
   - `length`：序列长度（正数）  

   示例：  
   ```r
   seq_data <- read.delim("seq_track.tsv", sep = "\t", stringsAsFactors = FALSE)
   ```  
   其中，`seq_track.tsv` 文件格式如下（示例）：  
   ```
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

2. **比对数据（ribbon_data）**  
   数据框，包含BLAST比对结果（可通过`outfmt6`或`outfmt7`格式转换），必须包含以下列：  
   - `qaccver`：查询序列ID（需存在于`seq_data$seq_id`）  
   - `saccver`：目标序列ID（需存在于`seq_data$seq_id`）  
   - `length`：比对长度  
   - `pident`：序列相似度（百分比）  
   - `qstart`/`qend`：查询序列上的比对起止位置  
   - `sstart`/`send`：目标序列上的比对起止位置  

   示例：通过BLAST结果文件生成：  
   ```r
   # 读取单条BLAST结果
   read_blast <- function(file) {
     df <- read.delim(file, sep = "\t", header = FALSE, comment.char = "#")
     colnames(df) <- c("qaccver","saccver","pident","length","mismatches","gapopen",
                      "qstart","qend","sstart","send","evalue","bitscore",
                      "qcovs","qlen","slen","sstrand","stitle")
     df
   }
   
   # 批量读取并合并多个BLAST结果
   blast_files <- list.files(pattern = "*.o7")  # 替换为你的BLAST文件路径
   all_blast <- do.call(rbind, lapply(blast_files, read_blast))
   ribbon_data <- subset(all_blast, length >= 100)  # 过滤短比对（可选）
   ```  

   你可以使用以下脚本对示例序列进行BLAST比对：  
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


### 基础使用示例  
```r
# 加载包
library(ggchord)
library(ggplot2)

# 调用ggchord生成多序列弦图
p <- ggchord(
  seq_data = seq_data,                # 序列信息
  ribbon_data = ribbon_data,          # 比对数据
  title = "多序列比对弦图",           # 标题
  seq_order = c("MT108731.1", "MT118296.1", "OQ646790.1", "OR222515.1"),  # 自定义序列顺序
  seq_gap = 0.03,                     # 序列间间隙比例
  seq_orientation = c(1, -1, 1, -1),  # 序列方向（1正向，-1反向）
  ribbon_color_scheme = "pident",     # 按相似度着色
  ribbon_alpha = 0.7,                 # 连接带透明度
  axis_gap = 0.1,                     # 坐标轴与序列的距离
  rotation = 45,                      # 整体旋转45度
  show_legend = TRUE                  # 显示图例
)

# 显示图形
print(p)

# 保存为高清图片
ggsave("multi_sequence_chord.png", plot = p, width = 10, height = 10, dpi = 300)
```  


## 参数详情  

| 参数类别       | 参数名                  | 类型          | 默认值                  | 描述                                                                 |
|----------------|-------------------------|---------------|-------------------------|----------------------------------------------------------------------|
| **核心数据**   | `seq_data`              | data.frame    | -                       | 包含`seq_id`（序列ID）和`length`（长度）的序列信息数据框             |
|                | `ribbon_data`           | data.frame    | -                       | 包含比对信息的数据框（需含`qaccver`、`saccver`等必要列）             |
| **基本样式**   | `title`                 | 字符串        | "Multi-sequence Chord Diagram" | 图形标题                                                             |
|                | `show_legend`           | 逻辑值        | TRUE                    | 是否显示图例                                                         |
|                | `rotation`              | 数值          | 45                      | 整体图形旋转角度（度，逆时针为正）                                   |
| **序列布局**   | `seq_order`             | 字符向量      | NULL                    | 序列展示顺序（默认按`seq_data`顺序）                                 |
|                | `seq_labels`            | 字符向量      | NULL                    | 序列标签（默认使用`seq_id`）                                         |
|                | `seq_orientation`       | 数值向量      | 1                       | 序列方向（1正向，-1反向），支持单值统一设置或按序列单独设置           |
|                | `seq_gap`               | 数值/向量     | 0.05                    | 序列间间隙比例（[0,0.5)），控制圆弧间距                               |
|                | `seq_radius`            | 数值/向量     | 1.0                     | 序列圆弧半径，支持单值统一或按序列设置                               |
|                | `seq_colors`            | 颜色向量      | NULL                    | 序列圆弧颜色（默认自动生成，基于RColorBrewer的Set1）                 |
| **连接带样式** | `ribbon_color_scheme`   | 字符串        | "single"                | 连接带配色方案：`single`（单一色）、`query`（按查询序列）、`pident`（按相似度渐变） |
|                | `ribbon_colors`         | 颜色向量      | NULL                    | 配色方案参数（单值/向量/渐变色阶）                                   |
|                | `ribbon_alpha`          | 数值          | 0.6                     | 连接带透明度（[0,1]）                                                |
|                | `ribbon_gap`            | 数值/向量     | 0.1                     | 连接带与序列圆弧的径向距离                                           |
|                | `ribbon_ctrl_point`     | 向量/列表     | NULL                    | 贝塞尔曲线控制点（默认圆心），用于调整连接带形状                     |
| **坐标轴设置** | `axis_gap`              | 数值/向量     | 0.05                    | 坐标轴与序列的径向距离（支持负值向内）                               |
|                | `axis_tick_major`       | 整数/向量     | 5                       | 主刻度数量                                                           |
|                | `axis_tick_major_length`| 数值/向量     | 0.02                    | 主刻度线长度比例                                                     |
|                | `axis_tick_minor`       | 整数/向量     | 4                       | 每两个主刻度间的次刻度数量                                           |
|                | `axis_tick_minor_length`| 数值/向量     | 0.01                    | 次刻度线长度比例                                                     |
|                | `axis_label_size`       | 数值/向量     | 3                       | 刻度标签文字大小                                                     |
|                | `axis_label_offset`     | 数值/向量     | 0                       | 刻度标签偏移量（正值向外，负值向内）                                 |
| **调试**       | `debug`                 | 逻辑值        | FALSE                   | 是否打印调试信息（有效/无效比对数量）                               |


## 图形解读  
- **序列圆弧**：每条彩色圆弧代表一条序列，长度按比例映射，箭头指示方向（正向/反向）。  
- **连接带**：连接不同序列的彩色区域，代表比对区间：  
  - 按相似度着色时，颜色渐变反映序列一致性（如从蓝色到红色表示相似度从低到高）。  
  - 按查询序列着色时，同一颜色的连接带来自同一条查询序列。  
- **坐标轴**：每条序列外侧的刻度线与数字，标注序列位置（单位与序列长度一致），方便定位比对区域。  


## 依赖包  
- [ggplot2](https://ggplot2.tidyverse.org/)：绘图核心引擎  
- [RColorBrewer](https://cran.r-project.org/package=RColorBrewer)：颜色方案支持  
- [grDevices](https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/00Index.html)：颜色处理工具  


## 许可证  
本项目采用MIT许可证，详情见[LICENSE](LICENSE)文件。  


## 版本更新记录  
### v0.1.0（最新）  
- 支持多序列展示（从双序列升级），通过`seq_data`和`ribbon_data`分离管理序列与比对数据。  
- 新增序列方向控制、自定义顺序、间隙与半径调整。  
- 实现可定制坐标轴（主/次刻度、标签位置）。  
- 连接带支持3种配色方案（单一色、按查询序列、按相似度渐变）。  
- 支持整体图形旋转与调试模式。  

### v0.0.2  
- 新增双序列的“弧线模式”与“直线模式”切换。  
- 支持弧线曲率调整与间隙控制。  

### v0.0.1  
- 初始版本，支持双序列BLAST比对弦图可视化。  


## 贡献与反馈  
欢迎提交issue报告bug或建议新功能，也可通过Pull Request参与开发。使用中遇到问题，可参考示例代码或开启调试模式（`debug=TRUE`）排查数据问题。
