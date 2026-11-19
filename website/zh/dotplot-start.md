数据来源、展示子集与重建规则见[旗舰数据来源与构建](flagship-data.md)。本页直接加载随包输入，集中设置绘图参数。

本页先用一对序列的真实 HSP 画出双轴线段，再扩展到四序列综合矩阵。最小输入是完整序列长度和 HSP 四端点；无需注释表。

{{code:courseblockam}}

<a id="dotplot"></a>

## 绘制HSP真实端点 {#dotplot}

随包 PAF 表保留一对序列的 11 条已观测 HSP。两条原始长度决定坐标轴；规范化后的反向 target 区间保留降序。


每条线从`(qstart, sstart)`连接到`(qend, send)`，无需`geom_seq()`。
颜色、透明度、线宽和线型复用`link_*`映射与标度。默认Identity范围为0–100，
刻度为0/50/80/90/95/100；显式用户标度优先。

{{code:dotplotview}}

{{figure:viewdotplot}}

保留反向 HSP 的降序端点。逐区间排序会改变原始线段方向。

`x_seqs`和`y_seqs`分别选择和排列accession。`gap`使用原生bp/aa单位表示显示留白，
不代表省略源序列。轴方向接受1/-1或完整命名向量。显示翻转改变斜率，但不改来源表。
每条序列从源位置1起算；拼接后的区段保留各自accession和局部刻度标签。

{{code:dotplotreverse}}

用`map_ggchord_position(..., axis = "x"/"y")`转换原生位置与拼接坐标，
`inverse = TRUE`用于反查，显示间隙返回NA。聚焦后的原始来源坐标用`source = TRUE`，
轴标签默认仍是窗口相对位置。未知方向会警告、省略并报告；方向元数据与端点冲突时报错。
同一显式HSP ID重复展示需要所有参与层指定`allow_repeated=TRUE`；无ID的重复或自比对
原行保持独立。没有cg/cs/CIGAR或其他精细轨迹时，只展示HSP线段，不推断逐碱基错配和缺口。

<a id="matrix"></a>

## 多序列综合矩阵 {#matrix}

以下直接使用随包数据 `flagship_dotplot_example`，计数与标签表已准备完成。完整保留四序列的 32 条 HSP，不补造自比对、对称记录或缺失的 HSP。266 条 CDS 属于共用来源表，但不交给 Dotplot 的线段层。空白区域表示本次输入没有关系记录，不能据此认定不存在同源性。

上方列标题标明 query 序列及完整的 1–length 范围，左侧行标题标明 subject 序列范围；下方和右侧只显示各自序列的局部坐标。沿一列和一行的交点即可定位一对序列。每个带边框的矩形保留两条序列真实长度对应的宽高，淡色背景只帮助分区，不编码相似性。矩形间的留白是显示间隙。

该例的每个无序序列组合只记录一个 query→subject 方向。空格表示没有输入记录；轴标题、分区和每格计数直接使用随包数据。


{{code:flagshipdot}}

{{manualfigure:flagdotplot}}

<a id="direction"></a>

## 按来源记录检查方向 {#direction}

四端点的相对方向来自源记录；反转显示轴会改变视觉斜率，但不改原始方向。可以按序列对统计已有 HSP 的方向：

```r
hsp_direction <- flagship_hsp |>
  mutate(relative_direction = if_else(
    (qend - qstart) * (send - sstart) < 0, "reverse", "forward"
  )) |>
  count(qaccver, saccver, relative_direction, name = "HSP")
hsp_direction
```

<a id="mapping"></a>

## 查询精确坐标 {#mapping}

下面将首条 HSP 的 query 起点换算到拼接 x 轴。反查时用 `inverse = TRUE`，并明确 `axis`；显示间隙返回 NA。该转换不改变四端点。

```r
map_ggchord_position(
  flag_dotplot, flagship_hsp$qaccver[1], flagship_hsp$qstart[1], axis = "x"
)
```

<a id="support"></a>

## 边界与下一步 {#support}

按有向序列对分面时，用 ggplot2 的行列分面分别指定两端 ID。需要逐位置 gap 或 mismatch 时，见[字符对齐](aligned-sequences.md)。
