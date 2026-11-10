# 统一关系模型与分阶段迁移提案

本文件记录首轮基线快照；后续实际实现见 [v0.17.0 开发记录](DEVELOPMENT-v0.17.0.md)。

适用基线：main `33236efa98f308895368367ac9f96ea7f7337566`，任务书 v10。
本文为首轮设计，除明确标为隔离原型的部分，均未接入公开 API。

## 数据模型与最小内部边界

用户仍输入普通 data.frame，额外列及 list-column 原样保留。优先给既有引擎增加桥接层，
不立即替换所有几何计算。以下是内部职责名，不是准备导出的函数：

| 内部操作 | 输入 → 输出 | 不变量 |
| --- | --- | --- |
| normalize_sources | 序列、注释、双端记录 → 源表及来源注册 | 不修改真实长度、坐标、链、pident、单位或文件行 |
| register_instances | 源序列 + 显示声明 → 唯一实例表 | accession 是实体键；instance_id 是显示键；不复制科学输入 |
| route_relations | 源关系 + 实例 + 图层选择 → 唯一配对表 | Coord 前确定路由；保留 query/subject 角色 |
| project_intervals | 已路由端点 + 窗口/方向/长度策略 → panel 局部几何 | 不在投影阶段重新过滤或复制 HSP |
| measure_layout | panel/guide/字体/设备 → 物理上下文、标签与范围 | estimated/measured、失败和字体回退分开记录 |
| compose_panels | 子图最终物理区域 + 配对 → 全局几何 | 先确定子图布局，再计算全局连接 |

源序列表至少含 `accver,length,sequence_unit`；注释保留 `start,end,strand,type` 及来源；
关系保留六端点字段、`source_hsp_id` 和已有证据列。注册表通过旁表保存来源，不向用户表
强塞保留字段；用户已有同名列必须在入口明确处理，不能覆盖。

来源身份应由输入对象注册而非每个图层临时 row_number() 产生。文件来源键可包含输入文件
哈希、原始记录行与 importer 信息；用户自备表可声明稳定行 ID；无声明时只能保证该对象
注册生命周期内的身份，不能声称不同进程/重排后的任意同值记录都能自动识别。
完全相同的两条源观测仍是两条记录，不能用数值去重。`source_row` 现有含义是解析后层输入行，
不改名冒充跨层持久 ID；新增旁表把它关联至 source ID。

实例表保留 `instance_id,accver,declared_order,panel_id`；关系路由保留
`source_hsp_id,displayed_relation_id,q_instance,s_instance,layer_id` 和重复显示状态。
窗口另保留 source_accver、origin、end、orientation、unit。panel 原生单位与全局 mm
分别存储，不能覆盖已有导出几何。新增字段采用按名消费，保留旧字段兼容。

## A–B–A：先确定关系，后换坐标

以下只展示未来用户语法，**当前包仍拒绝重复 seq_order，也不接受 between**：

```r
p <- ggchord(seq_data = seqs, ribbon_data = hits) +
  geom_seq(seq_order = c(A_high = "A", B = "B", A_low = "A")) +
  geom_link_ribbon(
    data = hits |> dplyr::filter(pident >= 90),
    between = c("A_high", "B")
  ) +
  geom_link_ribbon(
    data = hits |> dplyr::filter(pident < 90),
    between = c("B", "A_low")
  )
```

`prototype.R` 提供独立 `register_instances()` / `route_hsps()` 来检验这一数据流；
不是在源码增加同名函数，也没有伪装 `p + coord_*()` 已可运行。输出两条源记录各路由一次，
命名/数字实例选择等价，有歧义直接错误，B–A 通道不会交换原始 q/s。filter 的 NA 处理沿用
dplyr；90 为示例作者的阈值。源表保留，绘图实例不增加科学样本量。

同一源行跨层重复出现时，注册阶段报出重复展示候选，默认阻止无声明重复路由；如何公开
允许重复展示仍待 API 评审，不新增强制规则表。空层保留 schema 及声明，不参与 scale
训练；原型暂只验证空结果返回，不代表完整空层契约已经实现。

参数解析顺序：实例级 > 实体级 > 通用默认，与命名向量出现顺序无关。未命名等长向量
按声明顺序，未知/重复名字报错；实体名与显式实例名若冲突，建议拒绝歧义声明。
Feature/Track/Axis/Label 默认继承源实体，覆盖仅影响显示。需要修改的索引包括
`ggchord-build.R` 的 seq_order 检查、`layout-collinear-transform.R` 的唯一 accession、
`sequence_reference`、`layout-track.R`、feature/gene/label/axis 的 accession 查表、
`collinear-utils.R` 的相邻对选择、层输入注册和布局导出。

首批 Collinear/Chord 共用注册与路由；Radial 完成后接同一数据流。Circular/Dotplot
不适配重复实例，误用明确报错。Sankey/Fusion 复用内核但有独立配置，不承诺 A–B–A
语法在这两种布局一键互换。

## Feature 输入与兼容

推荐保持第三位置参数 `gene_data`，在原有形参之后、`...` 之前增加命名的 `feature_data`，
避免改变前三个位置参数。两个入口同时提供直接错误，提示显式合并表或提供图层 data。
正常化只在构造边界进行；旧 gene_data 使用每次构造至多一次的迁移 message，build/print
不重复提示。提示级别、频率和最终删除版本在公开实现前定稿；现在不承诺删除日期。

`geom_*(data=...)` 优先于构造数据；旧无 type gene_data 保持原来的 gene 表语义。
新 feature_data 中 geom_gene 默认仅选明确 CDS/gene，要求已知链；建议首批不自动纳入
tRNA/rRNA，保留在 generic feature，避免把 RNA 当蛋白编码箭头。该类型集合是推荐方案，
尚未作为公开默认固化。混合表缺 type/role 时给可操作诊断，要求显式子集，不猜测。

geom_feature 仍接受全表。建议保持层独立，不通过“图里碰巧还有 geom_gene”静默改变其数据；
同时绘制时使用显式、互斥的 filter，并在诊断指出跨层重复源 ID。`prototype-features.csv`
以合成 CDS/tRNA/rRNA/promoter/repeat 验证显式分区无丢失/无双绘。这只验证分区，尚未验证
新构造入口、迁移消息、aa、无 type 兼容和图层 override 的公开实现。
真实 strand、显示方向与内部辅助方向必须独立；已有 `.feature_biological_strand` 可复用，
不可因内部将未知方向转换为几何辅助 + 而改写科学链。

## Track：统计空间、数值空间、物理 band 分开

基线 bar 是 StatIdentity 数值柱。首批迁移推荐在同一候选阶段实现通用 geom_track、
geom_track_col 和 bar/count，再迁移已计算数值的示例；不能提前改官网调用或仅改帮助。
旧数值脚本应明确改用 col/identity，显式 bar/identity 作为经验证的通用组合允许保留。
bar 不根据“恰好存在 track_value”猜旧语义；对冲突输入在入口给迁移错误，防止悄悄变成计数。
历史版本脚本和冻结生成物保留历史事实，不批量更名。

`track-bar-migration.csv` 列出基线调用与引用，逐项区分实现、usage、运行示例、历史验收脚本。
README 通过共享示例生成，修改上游 snippet；不得漏掉间接来源。原型仅比较源位置加权 count
与 ggplot2 原生结果，并未为包内 bar 增加 Stat。

推荐流水线：source/window/instance/group 分组 → Stat → staged mapping/数值 scale 与
ggplot2 Position → track band 投影 → 曲线细分 → Geom。transform 与 Stat 的精确先后需逐
Stat 对照 ggplot2，特别是 transformed scale 下 bin/smooth，不能简单统一成一个未经验证顺序。
count 按相同 pos 聚合，weight 求和，不自动分箱；区间事件需要明确位置锚定，首批建议只收
pos 点事件。bin 必须另定义 boundary/closed/binwidth 与 focus 裁剪口径，未定义前不开放。

候选公开层分离 `track=position_track(...)` 与数据 `position="stack"/"identity"/...`；旧
`position=position_track(...)` 只作兼容的 band 声明，重复声明冲突报错。局部 x/y 只在通用
Track 标准化，不修改全局画布映射。支持 after_stat 需要真实 Stat 字段及数值验证；after_scale
继续使用 scoped aesthetic 改写，但不能把改写存在当作所有 staged mapping 都已支持。

| 阶段 | Geom | Stat | 数据 Position | 视图与前置条件 |
| --- | --- | --- | --- | --- |
| 基线 | bar/line/tile | identity | position 参数被 band 占用 | Chord/Circular/Collinear；已有行为，不是新语法兼容 |
| 首批目标 | col/bar/line/point/tile | identity/count | identity；col/bar 的 stack 经比对后加入 | 先 Collinear 数值等价，再验证两类曲线局部法向 |
| 第二批 | path/step/area/ribbon/segment/rect/linerange/errorbar/text | identity/bin/summary | stack/dodge 经分组与窗口验证 | 明确 xmin/xmax/ymin/ymax、tessellation 和裁剪 |
| 后续评估 | 上述已验证组合 | smooth/density/ecdf | jitter 等显式随机性策略 | 可追溯统计输出、单位、seed 与性能通过后开放 |
| 不自动支持 | boxplot/violin/raster/facet/全局 coord_flip | 任意第三方 Stat | 任意未验证 Position | 错误说明，不作盲目透传 |

同 track ID 的 layers 共用 band 与值域，不按图层数量分配多份 band；order 是物理远近，
添加顺序是 z-order。group、track group、track ID 不互相替代。Track group 首批不公开。
每项新组合都需 Geom×Stat×Position×Coord 加输入列条件的独立证据，目前均是目标矩阵。

## Dotplot、辐条与镜像布局

Dotplot 将一行 HSP 投为 `(qstart,sstart)→(qend,send)`。`prototype-dotplot.csv` 包含正、负斜率；
`real-paf-dotplot-endpoints.csv` 来自已有真实 PAF，并验证 importer 的 +1 和 reverse target
端点规则。缺失 cg/cs/CIGAR 不推断逐碱基 mismatch。多序列需轴集合/offset/方向与可逆来源。

`geom_link_line()` 当前要求 qpos/spos，连接两条骨架上的点，不能用它把一行 HSP 变成二维
斜线而丢掉另一对端点。推荐增加一个通用 alignment mark 薄层，共用 pair 校验与来源，
配合新双轴 Coord；首批使用普通 GeomSegment 作内部几何论证，无需新增 Stat。

Radial 固定 spokes，不再比较同心环方案。角度为 `theta_i=theta0+2*pi*(i-1)/N`，
源边界按 proportional/equal 策略映为 `r_inner+k_i*u`，再转 x/y。输入关系按环形邻接
筛选，N=2 只保留一份实际记录；N=1 拒绝多序列 synteny。原型输出 AB/BC/CD/DA 四条
与去掉 DA 的三条，未生成额外 HSP。PNG 只画中点连线检查拓扑，**没有实现 Ribbon
宽度、反向扭转、轨道/标签或通用 Coord**。正式原型下一步必须检查真实双端边界。

Sankey 左右各完整 A/B/C，内部注册 left/right 实例，默认过滤同一实体自连接，不按 ID 对
删除 A→B/B→A，也不自动镜像。建议以逐行忠实显示作为 MVP；规范化去重及显式镜像
留作后续选项评审。流带两端宽度来自各自区间，不宣称流量守恒。

纵向 Collinear 扩展现有 Coord，候选 seq_axis 与已有 orientation 区分：前者是全局主轴，
后者是逐序列源方向。建立局部切向/法向框架，骨架、轴、标签、控制点、ribbon 一起投影，
不整体旋转文字 grob。需覆盖 align/anchor/segments；当前仅有设计，尚无纵向 build。

## equal 长度与首末坐标轴

默认 proportional，用户 opt-in equal；统一布局策略管理 D_i，不能每个 Coord 各搞一套参数。
候选入口需结合没有 geom_seq 的 Dotplot 设计后定稿。边界 u∈[0,L] 投影 D*u/L，闭区间
[start,end] 用 [start-1,end]，中心用 (pos-.5)；`prototype-equal-boundaries.csv` 检验边界，
但现有代码某些位置用 pos/L，正式落地前必须逐几何桥接并记录 1 bp 差异，不静默重写默认。

导出 source_length/display_length/factor/mode/unit/instance；Chord equal 是角配额，半径
不同不保证物理弧长相等。equal 下共享绝对 bp/aa 比例尺失效，首批拒绝全局比例尺或改为
明确分序列局部标尺，图注/诊断必须可发现归一化，真实数据不变。

轴默认候选为实际显示区间首末源位置，不直接把 0 改字成 1。L=1 合并；focus/segment
用其源窗口端点；用户 breaks=NULL、显式 breaks、关闭轴、formatter 优先。
中间 pretty 候选严格取开区间，以真实字体盒及局部轴方向确定保留子集。
闭圆 1/L 贴近时用错位/分层或组合说明，保留两个 identity，不将 L 称为 0。
axis_selection 补充候选、选择、原因与 measured 状态；该求解器尚未实现。

## Ribbon、标签、图例与画布

填充优先级：显式固定值 > 显式 mapping 及用户 scale > 未指定时端点序列色渐变。
显式 pident 且无自定义标度才用固定 0–100、0/50/80/90/95/100 五段；其他数值不套用。
裸 scale 本身不声明数据变量，不能从“用户加过 scale”推断 pident。
端点渐变读取最终骨架颜色，沿路径 t 插值；评估自适应切片与设备原生渐变，原生线性填充
不能证明曲线路径渐变。边框只画一次；target 箭头标出完整目标区间肩部，尖点不冒充唯一
匹配碱基，query/subject 与生物学转移方向分开。

基线 pident breaks 少 95，默认 guide 为 colourbar，不可仅补 tick 后声称非等宽色阶已完成。
原生 ggplot2 的 guide_coloursteps(even.steps=FALSE) 是待验证组合之一；还要验证实际切割、
端点归属、图例 50/30/10/5/5%、90/95/100 文本避让。`guide-prototype.R` 只保存候选颜色取值
和要求的比例表；尚未检查实际 guide 段长，因此不称为目标图例已验收。

标签基线存在 restriction≥30 或混合 primer 的 dense 分支，以及 Inf 不隐藏冲突标签。
这些是源码机制，不等于图库混乱的已证根因。下一实验按任务书 A.6 保持数据/字符串/画布
单变量，记录 total/visible/hidden、引线交叉/穿字、文字框/feature/裁切与最小字号；对照
同尺寸同注释的 SnapGene 输出，未做视觉对照前不宣称达到其质量。先复用已有 dense
排序/装箱和测量，统一注册全部类别与可回溯候选，避免只改 make_catalog_plot 参数。

GuideMeasure 复用 diagnostics 的 gtable 宽高，但需要单项 key+label、title、role、锁定维度。
离散 nrow/ncol 候选联合 panel 可用宽高；连续/分段条只调方向、尺寸与文字，不拆行。
多 guide 按占用、间距和 role 顺序联合求解，显式 theme/guide/隐藏设置优先，无法容纳报告
原因与建议尺寸，不删除类别。隔离实验对 20 类在一个 8×6 PNG 设备上实测候选，不涉及
多 guide 或用户锁定约束，不证明已接入 ggsave。

fixed 保留基线绝对尺寸语义；scale 限幅调整默认角色尺寸；adaptive 重新测量/排版。
建议对象上的策略组件供 print/ggsave/viewer 共用，主题只供样式；公开名暂不定稿。
实际 panel mm、guide 占用与字体为依据，不能只看 DPI/短边。显式 size/theme 不被自动覆写。
缓存签名需增加最终 panel、字体/设备类型、主题/guide、策略、实例/路由与源对象修订。
现有 render-context 只含 status、units_per_mm、canvas 与设备编号的键需继续审计，
不能当作完整字体签名。保留有限迭代、原设备恢复、对象/缓存隔离和失败诊断。

## Fusion、语义放大与短序列研究

Fusion 是组合层，不是 Coord：子图独立构建 → panel 最终物理区域 → 源区间映射到局部
几何 → 全局 mm 变换 → ribbon/leader。对同一 accession 多 panel/实例要求显式选目标。
F0 先一个圆图与一个线性图、一条真实双端连接；F1 四方位；F2/3 再评估八方位、复杂
路由和排版。当前单 panel measurement 限制是前置风险，不能只拼图片称功能实现。

放大区域用 source accver/start/end；复用 focus 的裁剪/来源记录，跨原点要显式 wrap，
detail 独立重绘。stat_scope 全局/局部必须明确，装饰 leader 与真实 HSP ribbon 区分。
可先单 inset，不依赖全部复杂 Fusion 完成；现有 ggmagnify 是二维区域机制参考，不能直接
把弧形基因组区间塞入画布矩形。未复制代码，不引入依赖。

短序列两级数据：A 局部真实字符串+source position；B 已对齐字符+alignment column+
gap/source 双向映射。只有长度或 HSP 端点无法产生字符和 mismatch。首原型需有来源的
60 bp DNA/40 aa 片段及真实对齐证据，验证负链互补、CDS phase、gap 和尺寸下限。
本轮未生成字符图，真实 gapped 对齐样本尚需选择；不编造科研样本。

官方资料复核（2026-11-10；仅做能力参考，不复制资产）：

| 项目 | 已核查能力与借鉴边界 | 来源 |
| --- | --- | --- |
| ggmsa | DNA/RNA/蛋白 MSA 与 ggplot2 图形组合；优先评估桥接 | [Bioconductor](https://bioconductor.org/packages/release/bioc/html/ggmsa.html) |
| ggseqalign | DNA/AA/字符串比对静态呈现，兼容 Biostrings/ggplot2 | [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/ggseqalign.html) |
| ggseqlogo | aligned sequences/频率矩阵，bits/probability；不重复自研字形引擎 | [作者教程](https://omarwagih.github.io/ggseqlogo/) |
| Jalview | 本轮目标官方 annotation 页面访问超时；列注释/共识能力暂沿任务书作待复核参考 | [待复核页面](https://www.jalview.org/help/html/features/annotation.html) |
| ESPript | 预对齐蛋白+相似性/二级结构的出版输出，不接入在线平台 | [官网](https://espript.ibcp.fr/ESPript/ESPript/) |
| Geneious | 注释、翻译与差异高亮可借鉴；编辑和比对算法不纳入 | [官方能力页](https://www.geneious.com/features/sequence-alignment) |

ggplot2 官方 [bar/col](https://ggplot2.tidyverse.org/reference/geom_bar.html) 和
[coloursteps](https://ggplot2.tidyverse.org/reference/guide_coloursteps.html) 支持本文语义比较；
真实支持仍以本地 4.0.3 源码与原型结果为准，不承诺任意上游 Geom/Stat 兼容。
