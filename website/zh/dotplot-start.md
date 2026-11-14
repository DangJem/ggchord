## 多序列综合矩阵 {#matrix}

以下输入来自仓库原始 FASTA 和标准 BLAST 表。完整保留四序列的 32 条 HSP，不补造自比对、对称记录或缺失的 HSP。266 条 CDS 属于共用来源表，但不交给 Dotplot 的线段层。空白区域表示本次输入没有关系记录，不能据此认定不存在同源性。

{{code:flagshipdata,flagshipdot}}

{{manualfigure:flagdotplot}}

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

## 查询精确坐标 {#mapping}

下面将首条 HSP 的 query 起点换算到拼接 x 轴。反查时用 `inverse = TRUE`，并明确 `axis`；显示间隙返回 NA。该转换不改变四端点。

```r
map_ggchord_position(
  flag_dotplot, flagship_hsp$qaccver[1], flagship_hsp$qstart[1], axis = "x"
)
```

{{code:courseblockam}}

## 绘制HSP真实端点 {#dotplot}

有来源记录的EMALE PAF序列对包含11条已观测HSP，文件声明的序列长度提供两条轴。
`read_paf()`将零基半开区间转换为内部一基约定，并将负向比对的target端点转为降序。
原始文件及来源见[输入下载](formats.html)。以下准备代码需要可选依赖`dplyr`。

{{code:dotplotdata}}

每条线从`(qstart, sstart)`连接到`(qend, send)`，无需`geom_seq()`。
颜色、透明度、线宽和线型复用`link_*`映射与标度。默认Identity范围为0–100，
刻度为0/50/80/90/95/100；显式用户标度优先。

{{code:dotplotview}}

{{figure:viewdotplot}}

旧构造验证器可能对这些合法反向HSP的降序端点给出警告。应保留已规范化的端点顺序；逐区间升序排序会抹去显示倒位。

`x_seqs`和`y_seqs`分别选择和排列accession。`gap`使用原生bp/aa单位表示显示留白，
不代表省略源序列。轴方向接受1/-1或完整命名向量。显示翻转改变斜率，但不改来源表。
每条序列从源位置1起算；拼接后的区段保留各自accession和局部刻度标签。

{{code:dotplotreverse}}

用`map_ggchord_position(..., axis = "x"/"y")`转换原生位置与拼接坐标，
`inverse = TRUE`用于反查，显示间隙返回NA。聚焦后的原始来源坐标用`source = TRUE`，
轴标签默认仍是窗口相对位置。未知方向会警告、省略并报告；方向元数据与端点冲突时报错。
同一显式HSP ID重复展示需要所有参与层指定`allow_repeated=TRUE`；无ID的重复或自比对
原行保持独立。没有cg/cs/CIGAR或其他精细轨迹时，只展示HSP线段，不推断逐碱基错配和缺口。

## 边界与下一步 {#support}

不支持骨架 gene/feature/track/ribbon 层、facet 或作为 Fusion 端点。需要 gap/mismatch 时，提供真实字符对齐，见[字符工作流](aligned-sequences.html)。操作示例：选一条线，找回四端点，分别解释斜率和配色。
