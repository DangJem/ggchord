# v0.16.1 绘图职责与可发现性核查

状态：本版已确认取舍的实施与有限验证记录；不是 v1.0.0 的长期 API 冻结。
起点为 `e43fbee`，验证脚本为 `tools/validate-grammar.R`；详细源码参数基线仍使用
`audit-generated/api_arguments.csv`，不把旧版生成文件冒充本版运行证据。
本版当前签名/职责和支持矩阵分别为 `grammar-parameters-v0.16.1.csv`、
`grammar-compatibility-v0.16.1.csv`，输入散列见 `grammar-source-hashes-v0.16.1.csv`。

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
| gene fixed | anno | horizontal | label_* 几何参数 |
| feature fixed | label | feature | label_* 包装，复用 gene 标签内核但默认语义不同 |
| gene repel | anno | radial | 连线及遮盖策略为 label_segment_* |
| feature repel | label | feature | 内部优先和紧邻轨道；external=FALSE 有圆图约束 |

line 的默认 avoid=none，NULL gap 回退 .035；ribbon 的默认 avoid=uniform。
straight line 拒绝 control point 与 branch。ribbon 的 NULL 控制点在 chord 中为
历史原点，在 collinear 中从端点推导；F34 已决定，不再重新提默认变更。
line/ribbon 的端点参数现在统一为 link_gap/link_ctrl_point，保留各自默认和单位。
gene/feature 的宽度统一为 width，固定形状统一为 shape；角色美学仍明确使用
feature_shape。gene 调用 feature 构造器，再绑定定向基因的严格输入与 gene_fill 角色。
共享内核不等于删除快捷入口或让未知链绕过 gene 的输入约束。

## dots、staged mapping 与错误边界

- constructor 拒绝未知 dots；geom 把允许的固定样式交给 ggplot2 layer/geom，
  未知样式可能由上游警告，不能承诺全包统一停止。
- theme 显式签名只列专属设置；通用文字/标题/图例通过 dots 命名传递，并参与角色
  继承。已经删除的 underscore 参数先给迁移错误。五个预设默认和标准覆盖直接对照一致。
- reader、数据转换和导出不属于本版绘图 dots 统一工作。
- after_stat 只能读取 Stat 暴露的计算列；after_scale 读取已映射的角色美学。
  使用 `ggplot2::stage()` 合并同一个美学的开始/计算/标度后表达式。
- 几何角色必须在布局前求值；late stage 和命名空间形式的 after_stat/after_scale
  均拒绝。纯 start 的 stage 不要求后期计算，仍可在输入阶段求值。
- 显式固定样式覆盖同名映射；省略固定值才使用映射。构建时依据调用者是否提供
  参数解析默认，继承映射不受默认固定值压制。被覆盖的视觉映射不再训练默认 scale。
- retired 参数可能先触发 R 的部分匹配歧义；不通过恢复旧参数来规避此错误，
  调用者须使用完整的新参数名。此行为不等于静默接受未知参数。

## 三坐标兼容证据

完整源码组合范围见 `audit-generated/grammar-compatibility-v0.16.1.csv`，其证据等级仍为
源码合同，不能整体升级为运行通过。本版独立运行记录见
`audit-generated/grammar-validation-v0.16.1.csv`。
已运行 chord/circular/collinear 的两个独立 gene 层及交换层/scale 顺序，确认独立
data/mapping/Position 和用户颜色；chord/collinear 的 density 组合确认 after_stat
与 staged after_scale 可求值。不是全部参数、任意 stat、Position 或设备的保证。

circular 拒绝 link/ribbon/stat；collinear 拒绝 restriction/primer/centre-label；
geom_seq_scale 仅支持 collinear。源码没有 guard 的组合不自动认证为支持。

## F39/F40 证据与处置

F39：修复前四个位点的 circular 输入中 radial/elbow/straight 返回相同几何。
用户选择独立分派后，straight 为直接路径，elbow 保留一个折点，来源行一致；
默认 radial 与修复前源码的直接对照一致。冗余 trunk 别名在本版删除并指向 radial。
直线不承诺绕开骨架或其他连线；本版未承诺全域无碰撞。

F40：theme 的公开检查拒绝等长命名角度向量；帮助也只承诺单个方向/角度。
因此 helper 的位置分支优先没有形成已复现的公开错误，本版保留公开标量域，
不扩大 theme 的向量格式，也不把该内部候选写成已修复坐标缺陷。

## 去留结论

保留 gene/plasmid 语义快捷入口、五个主题预设、英美拼写别名和22个重导出；
重导出各有独立工具职责，并无同义自有参数重复。帮助按映射/标题/外观/保存/grid
用途列出来源；未据词法引用数删除用户便利能力。删除旧 offset、重复固定形状入口、
角色前缀的同义参数及 trunk 兼容取值，具体旧/新调用与 v0.16.1 删除版本见迁移表。
