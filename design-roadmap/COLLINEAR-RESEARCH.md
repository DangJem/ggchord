# v0.14.0 共线性调研：gggenomes 与相关 R 包

调研日期：2026-10-28。代码基线：`codex/v0.14.0`，`e6d99fb`，开发版本
`0.14.0.9000`。当前正式版仍为 v0.13.1。

**性质：以下第 1–7 节保留调研时的结论、问题与待确认建议。** 后续用户已授权修复
及功能安排，当前取舍见第 8 节与 [v0.14.0 实现记录](versions/v0.14.0.md)。
调研阶段不修改绘图
实现，不变更版本范围，不运行包级检测或发布流程。既有轻量验证通过，不代表下列
新输入场景已经覆盖或修复。

## 1. 总体判断

当前版本已建立基本线性坐标、同套 geom/scale、Position、focus 与布局导出的
闭环。与 gggenomes 的主要差距在比较基因组工作流：数据身份与组织、比较对象选择、
方向证据、锚点对齐、边界语义和可逆显示变换。继续增加普通箭头形状的收益较低。

建议先修正确性和科学表达，再评估工作流功能。现有独立图层、role-specific
scale/theme/guide、显式预处理与三种坐标共存的设计应保留，不能为模仿另一包而
把数据集合、Position lane 和物理 track 混成一个概念。

## 2. 已复现的当前问题

以下均为本地小型直接复现，没有执行 testthat 或包级检查。优先级是本次审查建议。

| 编号 | 优先级 | 复现及影响 | 建议处理 |
| --- | --- | --- | --- |
| C1 | 高 | 合法 `seq_gap = .01` 使相邻骨架 y 为 0、-0.02，但默认 point link 两端位于 -0.035、+0.015；ribbon 亦越过对方轨道。 | 给连接端点、feature 带和相邻骨架建立共同空间预算；自动间隙无解时适配或明确诊断。显式 gap 不得无提示缩小。 |
| C2 | 高 | `focus_ggchord_data()` 仅裁切/平移顶层 start/end，不同步转换 `segments`；后续 renderer 从 segments 重建旧坐标，注释跑出新序列。 | 统一处理所有受支持的坐标字段，裁切后重建派生 runs；不能支持的复合输入应显式拒绝，不能返回看似合法的结果。 |
| C3 | 高 | 真实不重叠的短 feature 被分配到相同 lane，但箭头显示工厂扩宽后实际相撞。 | lane 分配使用最终显示包络，或使用保留真实区间的短箭头策略；显示扩宽与生物坐标分别记录。 |
| C4 | 中 | `seq_order = c("A", "A", "B")` 被接受；导出表记录两个 A 的不同 y，但绘图按重复名称取到第一个 A，出现幽灵轨道。 | 拒绝重复 ID；是否允许 subset 维持现有选择语义，不把本问题顺带改成新的完整排列要求。 |

代码入口：

- C1：[`R/layout-collinear.R`](../R/layout-collinear.R) 的 `make_front()`、
  point-link `endpoint()` 与 `compute_collinear_layout()`。
- C2：[`R/focus.R`](../R/focus.R) 的 gene 裁切，结合
  [`R/common-features.R`](../R/common-features.R) 的 segments/runs 展开。
- C3：[`R/position-feature.R`](../R/position-feature.R) 的 lane 分配与
  [`R/helpers-feature-geometry.R`](../R/helpers-feature-geometry.R) 的最小显示区间。
- C4：[`R/ggchord-build.R`](../R/ggchord-build.R) 的 seq_order 校验与线性布局命名索引。

### C2 的具体输入

序列 A 长 1000；feature 为 400–700，其两个 segments 为 400–500 和 601–700。
截取 450–650 后，顶层 feature 为 1–201，新序列长 201，但 segments 仍为原坐标。
骨架的 x 范围为 0–2，feature 却约为 3.98–6.97。它会影响已有复合注释和 focus
组合，不能仅以“以后再支持复杂转录本”解释。

### C3 的具体输入及图义风险

长度 1000 的 A，上面两个正链 CDS 为 100–101 和 102–103，使用现有
`position_feature_stack()`。结果均为 lane 0，显示 x 范围分别约为
0.19225–0.20975 和 0.19625–0.21375，y 范围相同，实际图形重叠。

扩宽本身是现有可见性取舍，不能一概视为程序错误，但其幅度需要可见且可追踪：
在 1 Mb 序列上，500000–500001 的 2 bp feature，默认箭头绘制宽度换算为约
8.75 kb，当前导出缺少明确的扩宽标记。若用户按箭长比较基因大小，容易误读。
gggenes 提供了另一种取舍：短基因压缩头部到真实区间，不扩大区间。是否改变
ggchord 默认策略应先确定保真与可见性的优先关系，不在调研中替用户作决定。
[gggenes 箭头契约](https://wilkox.org/gggenes/reference/geom_gene_arrow.html)

## 3. 当前还应改进的科学表达和工程边界

### 3.1 截断必须可见，显示裁切不等于重新比对

现有 focus 默认 trim；源 accession、窗口起止和 source row 有记录，不能说原始
来源完全丢失。但截断端未明确作为图形和导出状态呈现，默认轴也只展示重定位坐标。

建议保留 source start/end、显示 start/end、左右截断标志及原坐标变换，并提供
断裂端标识、相对/源坐标轴的明确选择。对 ribbon 的裁切，目前按两端坐标线性插值，
不能把裁切后的子区间、按比例调整的 length 或沿用的 pident 解释为重新计算的比对
结果；有 gap/CIGAR 时尤其应区分显示近似与精确投影。

gggenomes 的官方边界教程明确讨论了“截断基因被读成短基因”和“整体区间相关不
保证任意子区间相关”的风险；其 drop/keep/trim 各有代价，不能仅因为另一包默认
drop 就静默更改 ggchord 的默认行为。
[边界处理](https://thackl.github.io/gggenomes/articles/marginal.html)，
[截断端装饰](https://thackl.github.io/gggenomes/reference/geom_seq_break.html)

### 3.2 标签避碰尚未等于完整引线避障

线性 shared solver 检查文本框、feature 和骨架，再重建引线；裁线主要面对本标签
层文字。跨层引线穿字、穿实体和引线之间相交缺少完整的联合检查。这是源码审查
确认的协议缺口；本次没有将所有碰撞类型都构造为已复现失败图。

建议把文字框、实体包络和引线段共同进入最终线性检查，优先保持局部短引线；
超出当前轨道可用空间时给出隐藏或空间不足原因，而不是把标签无限推向邻行。
不要求照搬圆形外部分配器，也不承诺任意非平面关系都能完全消交叉。

### 3.3 输出物理尺寸与可逆坐标变换

当前最长骨架归一化为宽 2，y 由 seq_gap 决定。原生 ggplot2 附加层的 x 不能直接
假定为 bp；focus 后的 source origin 又是另一层转换。建议导出稳定的 source
origin、方向、显示平移和比例，并评估按 accver 进行正反转换的辅助入口。

箭头头部、实体厚度、文字和 padding 还应区分基因组长度与毫米/字号。gggenes 的
along/away 抽象和真实 viewport 内延迟单位换算值得研究；不能简单为了共享代码
把三种 coord 都变成同一半径模型。内部重构和物理单位新 API 成本较高，应与
v0.16.0 输出质量设计衔接。
[gggenes 内部坐标设计](https://wilkox.org/gggenes/articles/coordinate-transformations.html)

### 3.4 现有优化器不等于线性自动定向/对齐

`optimize_ggchord_layout()` 已有排序和方向推荐，不能列为“完全没有”。但其评分
仍以环周分配后的端点交替与 chord span 为目标（`R/ribbon-layout.R`），不能把
所得分数宣称为线性轨道交叉数或锚点对齐质量。

若增强，需明确 coordinate-specific objective：邻接连接方向支持、轨道跨度、
相同轨道对中的线性交叉、锚点位移分别评价；输出评分依据和未决定项，保持用户
显式应用结果。不能让 build 悄悄改变序列方向或运行外部比对。

## 4. 调研时的功能候选（当时全部待确认）

| 候选 | 具体收益与边界 | 借鉴来源 | 建议归属 |
| --- | --- | --- | --- |
| 显式锚点对齐 | 按用户选定基因/区间的左端、中点或右端对齐，仅改变显示平移，导出锚点与变换。缺失、多拷贝锚点须明确处理。 | [gggenomes align](https://thackl.github.io/gggenomes/reference/align.html)、[gggenes 对齐](https://wilkox.org/gggenes/reference/make_alignment_dummies.html) | 共线性工作流首选候选；当前非目标，纳入 v0.14.0 前须确认范围。 |
| 有支持度的自动定向 | 基于指定 links 推荐 orientation，给出正/反支持、阈值和冲突；不声称生物学方向或全局最优。 | [gggenomes flip/sync](https://thackl.github.io/gggenomes/reference/flip.html) | 与锚定分开设计，待确认；先明确现有优化器目标。 |
| 显式比较对选择 | 邻接/指定 pairs/全部连接可选，报告保留与省略数；不自动丢弃跨轨道证据。先做小型过滤工作流也有价值。 | [gggenomes links](https://thackl.github.io/gggenomes/reference/geom_link.html) | 可优先评估的线性可读性增强。 |
| 稳定 feature ID 与命名数据集合 | 标签筛选或重排后仍能按身份继承 geometry；避免仅靠整表坐标相等匹配。集合名与物理 track 严格分开。 | [gggenomes 数据提取](https://thackl.github.io/gggenomes/reference/pull_track.html) | v0.15.0 的数据契约候选。 |
| genome → contig 两级组织 | 当前一条 accver 一行；增加一个 genome 的多个 contig 同排、保留 contig 身份与间隙，适合碎片化组装。不得拼接伪造连续序列。 | [gggenomes 输入](https://thackl.github.io/gggenomes/reference/gggenomes.html)、[序列布局](https://thackl.github.io/gggenomes/reference/layout_seqs.html) | 数据与布局共同设计，成本高，不建议塞入收尾补丁。 |
| 目标驱动 locus 生成 | 按 feature 条件找目标、按距离合并、非对称上下游扩展，先标出候选再 focus；复用已有多 locus 输入。 | [gggenomes focus/locate](https://thackl.github.io/gggenomes/reference/focus.html) | v0.15.0 数据准备候选。 |
| 同源簇和蛋白比对适配 | feature ID/cluster ID 转为显式 gene links；蛋白局部比对投影要处理负链、phase 与剪接，不能一律乘三。 | [gggenomes add_tracks](https://thackl.github.io/gggenomes/reference/add_tracks.html) | v0.15.0 可选适配器候选。 |
| 父子基因模型/domain 覆盖 | 子区间继承父基因翻转、裁切、轮廓和对象级 lane，并保留独立填色。现有 segments 并不等于完整 transcript/domain 模型。 | [gggenes subgene](https://wilkox.org/gggenes/reference/geom_subgene_arrow.html)、[gggenomes 基因模型](https://thackl.github.io/gggenomes/reference/geom_gene.html) | 先在 v0.15.0 评估 ID/Parent/phase，几何扩展另行确认。 |
| 公共比例尺与选择性序列轴 | bp/kb/Mb 格式化可先给示例；共享 bp 比例尺、只显示部分行的轴，可减少重复文字。保留用户 scale 优先。 | [genoPlotR 官方教程](https://genoplotr.r-forge.r-project.org/vignette.php) | 低成本示例优先；新 API 待确认。 |
| 有字号下限的内部文字适配 | 允许有界缩放、物理 padding 和换行，放不下再转局部行/隐藏；不以不可读小字“解决”碰撞。已有自动标签能力不重做。 | [gggenes 标签](https://wilkox.org/gggenes/reference/geom_gene_label.html) | 与 v0.16.0 物理输出联动。 |
| 定量 coverage/GC/variant 轨道 | 内容、数值 scale、带宽、y 轴和障碍协议分别建模；现有 lollipop 形状不等于定量轨道。GRanges 等用可选适配，避免强引入整套依赖。 | [trackViewer 官方教程](https://bioconductor.org/packages/release/bioc/vignettes/trackViewer/inst/doc/trackViewer.html) | 保持既定 v0.15.0 互操作、v0.16.0 定量轨道边界。 |

## 5. 不建议立即扩张的方向

- 多窗口同排、局部翻转：genoPlotR 已有分段 xlims/offsets 的表达，但 ggchord 需先
  定义分段变换、断裂处连接裁切与原坐标回查。不能用一个新的 orientation 参数
  草率覆盖断裂拓扑。[genoPlotR 官方教程](https://genoplotr.r-forge.r-project.org/vignette.php)
- 任意 facet、自由 track 容器、交互式浏览：继续保持既定非目标，先稳定现有共享
  grammar 的行为。
- 自动运行比对、隐式下载或强绑定第三方包：保留显式数据准备原则。
- 为“功能数量”新增重复箭头、第二套线性 geom 或品牌化函数名：现有形状与层
  已有基础，优先完善工作流及正确性。

## 6. 建议执行顺序与验收意图

1. **当前 v0.14.0 修复候选**：C1–C4；优先确保数据坐标、可见几何与导出一致。
   只为确认缺陷保留最小回归，复杂质量矩阵留在显式工具脚本。
2. **科学表达收敛候选**：截断/扩宽标记、原坐标来源、引线最终检查、丢弃/隐藏原因。
   本次仅记录；默认短箭头与边界策略若需改变，另行确认并记录迁移。
3. **功能取舍待确认**：先评估锚点对齐、方向支持和比较对选择。它们应共享一个
   可逆布局变换表，不宜分别新增相互独立的“调图参数”。
4. **后续版本候选**：稳定 ID、命名集合、格式及蛋白/同源簇适配归数据准备设计；
   物理尺寸和定量轨道归输出/轨道设计。不因本次调研自动调整 v0.15/v0.16 范围。

验收关注公开语义：改变显示方向后所有相关层一致；过滤或重排标签不改变实体；
裁切后所有顶层及嵌套坐标一致；图形是否被放大/截断可追踪；自动推荐有证据与
不确定性；显式选择/隐藏均可报告。不得以精确内部浮点坐标、迭代次数或固定秒数
冻结实现。

## 7. 调研术语与来源边界

- `sync()` 是基于连接支持的自动方向翻转；`align()` 才是目标驱动水平锚定。
  既有规划中的“类似 sync 的自动对齐”应按此区分；本报告不借术语修正扩张版本承诺。
- gggenomes 的主参考页标注 1.2.0，个别已读取参考页仍标注 1.1.3；gggenes 的已读
  页面标注 0.7.0。按各页具体接口讨论，不把所有能力宣称为新版本新增，也不据此
  推断 CRAN 最新发布版本。
- 来源为维护者官方文档、官方仓库和 Bioconductor 官方文档；仅研究设计与公开契约，
  未复制第三方实现、资产，也未将任何新依赖加入 ggchord。


## 8. 用户授权后的执行安排与当前状态

用户随后授权自主修正与增强。以下记录取代本报告中相应的“待确认/当前非目标”
判断；第 1–7 节仍保留基线 `e6d99fb` 的发现，不能据此推断缺陷仍存在。

| 调研事项 | 已确定安排 |
| --- | --- |
| C1–C4 | v0.14.0 已实现通道预算、focus 复合坐标转换、保留真实区间的短箭头及重复顺序校验。 |
| 锚点、方向证据、比较对 | v0.14.0 已实现 `coord_collinear(align, anchor, orientation, min_support)` 及 ribbon 邻接/pairs 选择；全部使用普通 `+` 组合，不增加 sync 变换步骤。 |
| 截断与可逆显示变换 | v0.14.0 已实现源坐标和截断记录、源坐标轴、`geom_seq_break()` 及 `map_ggchord_position()`；显示裁切不声称重新比对。 |
| stable feature ID、引线 | v0.14.0 已实现 ID 继承筛选/重排标签的 lane、联合引线检查和隐藏原因；命名数据集合安排在 v0.15.0。 |
| 生态互操作 | issue #2 的 Seqinfo/GRanges/GInteractions 最小可选适配提前纳入 v0.14.0，含通用转换器和适合的直接输入；更多格式/生态工作流仍在 v0.15.0。 |
| genome/contig、目标 locus、同源簇/蛋白、父子模型 | 已纳入 v0.15.0 数据契约和适配安排。 |
| 多窗口同排、局部翻转 | 已安排在 v0.15.0 分段变换设计中，保留断裂和源坐标语义后实施。 |
| 物理尺寸、最小字号、公共比例尺、选择性轴、定量轨道 | 已纳入 v0.16.0 输出与轨道设计。 |

本轮是开发实现与有针对性的轻量验证，未正式发布；具体实现、验证范围与待发布
验收见 [v0.14.0 版本文档](versions/v0.14.0.md)。相关设计依据仍来自第 1–7 节官方
来源；实现未复制第三方代码或资产。Bioconductor 包现在仅列为可选 Suggests。
