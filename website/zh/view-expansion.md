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

## 真实比例的径向辐条 {#radial}

随包四序列数据保留原始长度和CDS链。以下显式展示筛选选择长度至少900 bp、identity
至少75%的已观测HSP，并取每个有向序列对最长两条。计数轨道统计所选query起点观测数，
不将其解释为深度或覆盖度。

{{code:radialdata}}

{{code:radialview}}

{{figure:viewradial}}

辐条顺序来自`geom_seq(seq_order=...)`。默认按真实长度比例展示；只有显式
`length_mode="equal"`才归一显示长度，源坐标保持不变。`start_angle`从正x轴逆时针计角，
`clockwise`控制辐条顺序，不改生物链。特征和轨道的正偏移指向辐条逆时针法向。
物理Track宽度和源统计复用共享band与Stat系统。

首尾序列只有在**存在HSP**时才在闭合通道连接；没有该关系时通道保持空白。
N=2每条原始HSP只绘一次，N=1报错。反向HSP保留完整、有顺序的前缘及扭转。
新Radial采用闭区间边界`[min-1,max]`，所以单位置HSP仍有宽度；轴位置仍是一基。
已有视图继续采用各自文档约定。中心孔只是显示留白，不删除序列；序列名位于辐条外端。

{{code:radialrotate}}

Ribbon可用显式`pairs`或`adjacent_only=FALSE`选择非邻接关系；非邻接点连接须指定
`between`。命名A–B–A实例沿用[明确路由契约](foundation.html#instances)。
用户控制点使用最终x/y布局单位。端点渐变和subject尖头仍须显式选择，默认保留Identity填充。

## 导出与反向映射 {#export}

{{code:radialmap}}

导出保留`radial_frames`、原始源区间与`radial_endpoints`中的完整投影前缘，以及连接
筛选、实例路由、轨道统计和标签可见性。Radial反查必须给出所选骨架上的x、y，避免
仅凭多条辐条共享的单个x值猜位置。focus/segment输入保留源原点；径向显示不将多个
窗口拼为一条辐条。

## 支持组合 {#support}

| 图层或操作 | Dotplot | Radial |
| --- | --- | --- |
| `geom_alignment_segment()` | identity HSP线段 | 报错 |
| 序列、gene/feature、region、break | 报错 | 源空间局部投影 |
| 手动/自动gene与feature标签 | 报错 | 共享实测队列；按目标尺寸检查隐藏 |
| `geom_track()`及已有Track快捷层 | 报错 | 共享源Stat与band投影 |
| Ribbon和点连接 | 报错 | 默认环形邻接或显式路由 |
| 原生ggplot2注释 | 全局拼接单位 | 最终布局x/y单位 |
| 共享branch、facet | 报错 | 报错 |
| 圆图专有引物、酶切、中心标签 | 报错 | 报错 |
| 共享直线基因组比例尺 | 报错 | 报错；使用每条辐条的原生轴 |

这些是已支持契约和有限示例，不代表穷举全部参数组合。中心通道无法容纳间距时会给出
可操作错误；显式控制点可能穿过中心或其他辐条。密集标签可按图层策略隐藏，导出保留
来源身份。复杂Fusion和精细比对轨迹属于后续阶段。
