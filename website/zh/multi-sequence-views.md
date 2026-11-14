## 一次准备 {#input}

相同来源输入不代表默认展示完全相同的关系。本例先明确全对选择再比较几何；随后核对各视图默认筛选。依赖：ggchord、dplyr。

{{prepare:radialdata}}

{{code:coursecomparison}}

### Chord {#chord}

![chord：同源数据全对关系](../assets/coursechord.png)

### Linear {#linear}

![linear：同源数据全对关系](../assets/courselinear.png)

### Radial {#radial}

![radial：同源数据全对关系](../assets/courseradial.png)

### Sankey {#sankey}

![sankey：同源数据全对关系](../assets/coursesankey.png)

## 默认选择与可比条件 {#selection}

| 视图 | 默认关系选择 | 本例条件 |
| --- | --- | --- |
| Chord | 全对 | 同一 HSP 表，显式全对 |
| Collinear | 相邻行 | 显式 adjacent_only = FALSE |
| Radial | 环形邻居 | 显式 adjacent_only = FALSE |
| Sankey | 排除同实体自身关系 | 本例真实表不含自身关系 |

{{code:coursecomparisonreport}}

## 保留旗舰与派生视图 {#flagship}

弦图和共线性图继续使用现有旗舰版本，保留构图、输入筛选、配色和图层设计；图像仅随本次已确认的默认线型缺陷修复而重绘。Radial 从现有弦图派生，Sankey 从现有共线性图派生。两套既有旗舰共用来源准备，但各自保留原有特征与跨行 HSP 的显示选择，因此不把四图说成完全相同的显示子集。

{{prepare:sequencepalette,docstyle,showcasedata,chordshowcasegenes,showcasechord,showcasecollinear}}

Radial 移除 Chord 专属的半径和曲率；Sankey 使用各列原生坐标轴，不照搬 Collinear 的 Cartesian 控制点和共享比例尺。完整图保留所有原 CDS 名称候选、序列图例和来源记录。注释在最终投影坐标中按真实文字尺寸避让，优先邻近特征的短引导线；数据密集时仍可能隐藏无法安放的名称，可用 `diagnose_ggchord()` 查看实际导出尺寸下的显示与碰撞结果。

{{code:flagshipmulti}}

### Chord {#flagship-chord}

{{manualfigure:flagchord}}

### Collinear {#flagship-collinear}

{{manualfigure:flaglinear}}

### Radial {#flagship-radial}

采用 12 × 10.5 英寸画布，辐条从 45° 起排、序列名水平显示。保留全部 20 条原 CDS 几何与名称候选，名称随源特征就近安排，避免全局两侧名列和长引导线压缩序列主体。Sequence、Strand 与按数值区间比例显示的 Identity 图例集中在下方。

```r
view_ggchord(flag_radial, width = 12, height = 10.5, units = "in")
```

{{manualfigure:flagradial}}

### Sankey {#flagship-sankey}

四条序列在两列中纵向排列，使用 10 × 16 英寸纵向画布。16 条原 CDS 几何与名称候选、10 条已选择 HSP 均保留；左右副本分别标注，名称不再每序列只选一条。列外的短连接与局部避让保持文字靠近对应特征；完整刻度和三组图例保留。首页预览另省略密集名称和坐标文字。

```r
view_ggchord(flag_sankey, width = 10, height = 16, units = "in")
```

{{manualfigure:flagsankey}}

## 方向、实例与显示比例 {#direction-and-adjacency}

生物学链、HSP 端点方向、显示 orientation 和 seq_order 各有独立含义。Sankey 的左右副本不增加观测，区间宽度不是流量守恒。图层仅在共同支持范围内可复用。

## 显示长度与连接样式 {#projection}

`length_mode="equal"` 改变角度分配或直线跨度，保留原始长度、区间与局部轴。
`sequence_display` 导出源长度和投影系数。多序列等长显示不能使用共享绝对比例尺。
默认轴显示 1 与真实末位，中间刻度按可用物理空间选择；闭圆接缝合并两端。
显式 breaks 和 labels 优先。

{{prepare:foundationdata}}

{{code:foundationequal}}

{{figure:foundationequal}}

默认 Ribbon 使用 Identity 填充并保留完整区间端面；显式 `geom_link_ribbon(gradient = TRUE)` 才沿曲线插值端点序列色。`tip="subject"` 才在 subject 端装饰性收尖，
箭肩保留平移后的目标端点跨度，尖点不代表唯一匹配碱基或生物学转移方向；
源端点仍在布局导出中。`tip="none"` 保留完整端面。固定 fill、ribbon_fill 映射或显式
Ribbon 填充标度会覆盖显式启用的端点插值。显式 pident 使用固定 0–100 域，边界为
0/50/80/90/95/100，五段数值比例为 50/30/10/5/5%。用户标度与 guide 优先。

### 可选 subject 收尖 {#tip}

{{code:foundationarrow}}

{{figure:foundationarrow}}

自动离散图例测量文字后选择行列，渲染循环联合测量图例与 panel；显式 guide/theme
设置优先。圆图外侧标签在所有数量下共享测量队列，`label_visibility` 报告源行、可见性
及可用的放置字段。未测得碰撞数不等于零。过小画布仍可能拥挤，请按目标导出尺寸
检查[输出诊断](output.html)。

自动 Identity 标尺按边界文字的实际尺寸预留长度，不再把通用色块尺寸当成整条标尺长度。手动指定长度时，使用 `legend.ribbon.key.height`、`legend.ribbon.key.width` 或显式 guide。圆图 parallel 坐标轴文字逐字符沿圆弧排布；水平和数值角度仍保留直线文字。

## 操作示例 {#exercise}

选一条 HSP，指出两端来源；恢复默认邻接后说明为何可能省略。按目标尺寸保存选定视图，保留输入表。
