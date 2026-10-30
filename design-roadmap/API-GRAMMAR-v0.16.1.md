# v0.16.1 绘图职责与可发现性核查

状态：开发中的证据记录；公开取舍尚待用户选择，不是已冻结的新契约。
起点为 `e43fbee`，验证脚本为 `tools/validate-grammar.R`；详细源码参数基线仍使用
`audit-generated/api_arguments.csv`，不把旧版生成文件冒充本版运行证据。

## 职责与单位

| 任务 | 入口/拥有者 | 单位与边界 |
| --- | --- | --- |
| 提供序列与默认业务表 | `ggchord()` | length/start/end 为 bp；不把窗口空隙补为原始序列 |
| 数据角色与视觉变量 | `aes()` | accver/start/end/strand 等角色在输入表求值；视觉映射在展开几何求值 |
| 值域、颜色、刻度与 guide | `scale_*()`、`guide_*()` | 用户 scale 优先；breaks 是原始 bp，不预先取 log |
| 基因/通用注释实体 | `geom_gene()`、`geom_feature()` | 生物学链与显示 orientation 分开；未知链不自动变正链 |
| 局部空间位置 | `position_strand()`、stack/plasmid/track Position | offset/width 使用局部布局距离；lane 是编号，非物理单位 |
| 点连接与区间连接 | `geom_link_line()`、`geom_link_ribbon()` | qpos/spos 与区间端点均为 bp；gap 是布局距离，控制点是布局坐标 |
| 连接统计 | `stat_ribbon_bundle()`、`stat_ribbon_density()` | 预先计算再经 Stat 暴露 bundle_n/bundle_weight/density；聚合来源不等同单一输入行 |
| 画布模式 | `coord_chord()`、`coord_circular()`、`coord_collinear()` | ratio 是 y/x 比；limits 是布局坐标；rotation/gap（圆图开口）为度 |
| 非数据外观 | `theme_ggchord()` 及预设 | 轴间距为绝对 grid unit；文字方向公开接受单个关键词或有限角度 |
| 独立 bp 参照条 | `geom_seq_scale()` | 仅 collinear，长度为 bp；offset/tick_height 可为布局距离或绝对 unit |
| 序列位置刻度规则 | `scale_seq_position_continuous()` | 三坐标序列轴的 breaks/labels/transform；不会新增独立比例尺条 |

`geom_seq(seq_gap=...)` 的单位随模式有明确差异：chord 为周长比例，collinear 为
最长显示骨架宽度的比例；circular 的开口由 `coord_circular(gap=...)` 以度控制。
collinear 拒绝显式 radius/curvature，不把环形参数静默解释成直线距离。
基因箭头 head length/width 与设备毫米线宽不是同一种尺寸；不把所有几何参数
统一为物理单位。本版不扩展 aesthetic、scale、guide 或 theme 的列表语法。

## 标签与连接现状对照

| 层 | 文本角色 | 默认布局/方向 | 说明 |
| --- | --- | --- | --- |
| gene fixed | anno | horizontal | gene_label_* 几何参数 |
| feature fixed | label | feature | label_* 包装，复用 gene 标签内核但默认语义不同 |
| gene repel | anno | radial | 连线及遮盖策略为 gene_label_segment_* |
| feature repel | label | feature | 内部优先和紧邻轨道；external=FALSE 有圆图约束 |

line 的默认 avoid=none，NULL gap 回退 .035；ribbon 的默认 avoid=uniform。
straight line 拒绝 control point 与 branch。ribbon 的 NULL 控制点在 chord 中为
历史原点，在 collinear 中从端点推导；F34 已决定，不再重新提默认变更。

## dots、staged mapping 与错误边界

- constructor 拒绝未知 dots；geom 把允许的固定样式交给 ggplot2 layer/geom，
  未知样式可能由上游警告，不能承诺全包统一停止。
- theme 的 dots 只进入 ggplot2 theme；已经删除的 underscore 参数先给迁移错误。
- reader、数据转换和导出不属于本版绘图 dots 统一工作。
- after_stat 只能读取 Stat 暴露的计算列；after_scale 读取已映射的角色美学。
  使用 `ggplot2::stage()` 合并同一个美学的开始/计算/标度后表达式。
- retired 参数可能先触发 R 的部分匹配歧义；不通过恢复旧参数来规避此错误，
  调用者须使用完整的新参数名。此行为不等于静默接受未知参数。

## 三坐标兼容证据

完整源码组合范围见 `audit-generated/api_compatibility_matrix.csv`，其证据等级仍为
源码合同，不能整体升级为运行通过。本版独立运行记录见
`audit-generated/grammar-validation-v0.16.1.csv`。
已运行 chord/circular/collinear 的两个独立 gene 层及交换层/scale 顺序，确认独立
data/mapping/Position 和用户颜色；chord/collinear 的 density 组合确认 after_stat
与 staged after_scale 可求值。不是全部参数、任意 stat、Position 或设备的保证。

circular 拒绝 link/ribbon/stat；collinear 拒绝 restriction/primer/centre-label；
geom_seq_scale 仅支持 collinear。源码没有 guard 的组合不自动认证为支持。

## F39/F40 证据与待定处置

F39：四个位点的 circular 输入中 radial/elbow/straight 返回完全相同的非空
restriction_sites 几何，符合源码中未读取 params$leader 的事实。兼容别名澄清与
实现独立路由两方案待用户选择；不提前修改默认或承诺实现不同线路。

F40：theme 的公开检查拒绝等长命名角度向量；帮助也只承诺单个方向/角度。
因此 helper 的位置分支优先没有形成已复现的公开错误，本版保留公开标量域，
不扩大 theme 的向量格式，也不把该内部候选写成已修复坐标缺陷。
