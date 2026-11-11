## 来源表与显示副本 {#source}

本章使用随包序列长度、真实观测 HSP 和源链注释。`feature_data` 是统一注释来源，
通用特征层继承全表；默认 gene 层要求明确的 gene/CDS 类别。图层显式提供的表或
数据函数优先。组合两类几何时请使用互斥子集。

{{code:foundationdata}}

prepare、clean、focus 和 segment helper 也接受并返回 `feature_data`，不能同时提供
`gene_data`。未知链仍为未知；旧 gene 输入保留方向语义和位置参数。

## 命名 A–B–A 路由 {#instances}

三个显示行来自两个生物学序列。`seq_order` 的值标识来源，名称标识实例。
先筛选每层的原始 HSP，再选择路由；`between` 接受两个名称或两个声明序号。
歧义路由会报错。以下两个子集互斥，合起来保留所有 HSP。

{{code:foundationinstances}}

{{figure:foundationvertical}}

同一个 `foundation_plot` 可以使用 `coord_chord()` 或 `coord_collinear()`。
纵向 Collinear 改变全局视图，保留源链和源坐标，局部正法向朝左。显式文字样式优先。
`map_ggchord_position()` 接受实例名称，纵向反查传入最终 y。

`export_ggchord_layout()` 提供 `display_instances`、`instance_routing` 和
`instance_omissions`。显式 `source_hsp_id` 可检测跨层重复显示；只有确实需要重复时才
设置 `allow_repeated=TRUE`。没有 ID 时身份只在层内有效，不能以行值相同认定同一来源。

## 镜像双列关系 {#sankey}

Sankey 将每个来源在左右两列各显示一次。左侧 query 端点连接右侧 subject 端点，
互反输入行保持独立；默认排除同实体自身连接，`self_links=TRUE` 可启用。
宽度由两端真实区间决定，表达的是关系，不是流量守恒。本版由坐标系管理左右副本，
不接受 `between`，也不支持共享分支。

{{code:foundationsankey}}

{{figure:foundationsankey}}

## 统计观测数或显示给定数值 {#tracks}

`bar` 按源位置统计观测数，`col` 显示已提供的数值。各 accession/window 在视图投影或
显示复制前独立计算 Stat 和数据 Position。`track=position_track(...)` 声明物理 band，
`position` 独立选择 stack、dodge 等原生数据调整。旧的
`position=position_track(...)` 物理 band 写法继续接受。

{{code:foundationcount}}

{{figure:foundationcount}}

下面的快捷入口在相同物理轨道中统计同一组观测。

```r
p_count_shortcut <- ggchord(foundation_seqs) +
  geom_seq() +
  geom_track_bar(
    aes(x = pos), data = foundation_events,
    track = position_track("query_starts")
  ) +
  coord_collinear()
```


| 几何 | 默认 Stat | 输入含义 |
| --- | --- | --- |
| bar | count | x/pos 观测，可选 weight |
| col | identity | 已给定 y/track_value，区间或位置 |
| point、line、path、step、text | identity | 源 x 与局部 y，text 还需 label |
| tile、rect、area、ribbon、segment、linerange、errorbar | identity | 对应 ggplot2 原生范围字段 |

可用统计包括 identity/count/bin/summary，bin 必须显式给出 binwidth 或 breaks。
支持原生 group 与 after_stat 映射；数值轴保留统计量单位。任意第三方 Geom/Stat 插件
及 smooth/density/ecdf 不属于首批支持矩阵。既有 track line/tile 快捷层保留其原有行为。

源位置使用 x orientation，原生文字字号目前使用 mm，区间几何只支持完整轮廓。其他方向、字号单位或局部轮廓会明确报错。原生文字样式和路径箭头继续转发。

首批原生 Track 要求单面板和 ggchord 坐标系。分面、`coord_flip()`、raster、boxplot 和 violin 尚不支持，显式请求会报错。纵向投影使用 `coord_collinear(seq_axis = "vertical")`。统计结果在 dodge 等数据位置调整前保留源成员。矩形和 tile 使用输入边界，不虚构零基线。

## 显示长度与连接样式 {#projection}

`length_mode="equal"` 改变角度分配或直线跨度，保留原始长度、区间与局部轴。
`sequence_display` 导出源长度和投影系数。多序列等长显示不能使用共享绝对比例尺。
默认轴显示 1 与真实末位，中间刻度按可用物理空间选择；闭圆接缝合并两端。
显式 breaks 和 labels 优先。

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
