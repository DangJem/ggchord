# API 目标架构与边界

**状态：v0.16.0 审计基线、v0.16.1 绘图语法及 v0.16.2 数据契约实施边界；长期冻结仍待后续阶段。**

沿用 [`AGENTS.md`](../AGENTS.md) 的长期规则。本文件记录当前源码观察、目标边界
及待决策设计；具体问题和证据见 [决策台账](API-AUDIT-DECISIONS.md)。

## 现有架构观察

- `R/ggchord.R::ggchord()` 承载 seq/ribbon/gene 三表、验证结果和 plot-owned ref，
  构造默认坐标与主题；`+.ggchord` 对新组合建立独立 ref，并处理用户 scale 优先。
- `R/zzz.R::ggchord_capture_layer_input()` 保存每层原数据和映射，
  `ggchord_resolve_layer_input()` 处理回退输入。构造器并不是所有图层的数据入口。
- `R/ggchord-build.R::compute_chord_geometry()` 按 layer_id 分组，构建 layer_inputs /
  layer_geometry，并组织障碍物、轨道和跨图层标签；不是只有一个最终图层状态。
- `R/layout.R::compute_chord_layout()` 仍在一个函数环境内依次 eval 各 quote 阶段。
  这是需要显式依赖图的现状，单凭这一结构不能认定几何结果错误。
- Coord 类提供不同坐标行为；`ggplot_build.ggchord` 和
  `ggplot_gtable.ggchord_built` 将布局与设备测量接入普通 ggplot2 流程。
- `R/render-context.R` 有模块级临时上下文 `.ggchord_render_state`，由
  `ggchord_with_render_context()` 通过 on.exit 还原；另有 plot-owned 缓存。
  临时动态上下文不是持久缓存，是否存在异常或嵌套泄漏仍需定向验证。
- `get_chord_layout()` 返回布局列表，`export_ggchord_layout()` 整理成分层表格。
  后者公开帮助已经描述 annotation_registry 等字段，不能未经迁移将其降为内部。

以上是静态源码观察；第一轮运行验证只覆盖迁移文档列出的修复和代表性圆图 build，
不等同于三坐标、所有图层或所有设备的组合验收。

## 目标数据与生物学契约

| 数据角色 | 目标边界 | 不应混淆的对象 |
| --- | --- | --- |
| sequence | 原样保留 accver、明确全长和坐标单位 | contig、显示窗口与物理轨道 |
| feature/gene | 区间、链、分段和原注释来源可追溯 | 生物学链与显示方向；gene 与任意 feature |
| link | 原样 query/subject、两端来源与关系 | 通用 link 与带特定统计字段的 alignment |
| quantitative observation | 数值、采样位置、缺测及来源 | 数值 scale 与 band 空间位置/设备宽度 |

所有转换都应明确 1-based inclusive 与外部 0-based half-open 的关系。跨原点、
单碱基、零宽、反向端点、越界和分段必须分别定义，不能通过猜测补造业务数据。
`strand` 表示原始链，`orientation` 表示显示选择，连接相对方向按各 API 明确计算；
不能只统一单词却改变既有数据含义。

read、as、validate、clean、prepare、focus/segment/project 的职责分别是读取、
适配、检查、显式修复、编排、选择/分段/投影。prepare 可以调用这些步骤，但必须
公开步骤顺序和报告，不隐藏不可逆变换。F03 在 v0.16.2 保留任务相关默认，方向敏感管线使用 keep；当前逐入口合同见
[数据契约](API-DATA-CONTRACTS-v0.16.2.md)。

## 绘图职责

v0.16.1 的职责、单位与实际运行边界见 [绘图语法核查](API-GRAMMAR-v0.16.1.md)。
本版已直接验证三坐标的独立 gene 双层、feature 与固定/repel 标签组合，以及
chord/collinear 的 density staged mapping；这不升级其他源码矩阵项的证据等级。

| 层级 | 目标职责 | 审计重点 |
| --- | --- | --- |
| constructor | 基础数据与全局验证 | 是否需要新容器仍待真实用例，不先加入口 |
| aes / scale / guide | 映射、值域规则与图例 | 固定值与映射分离，用户 scale 始终优先 |
| geom | 可识别的视觉对象与固定外观 | 保留每层 data/mapping/参数，便利封装可共享内核 |
| stat | 有明确含义的计算与聚合 | after_stat、来源集合与数据损失说明 |
| position | 序列局部空间位置 | 不能承担全画布坐标或数值 scale |
| coord | 全画布模式及坐标转换 | chord/circular/collinear 差异显式，不伪造兼容 |
| theme | 非数据元素及外观预设 | 基础元素与真正自定义元素分别说明 |
| preview / export / diagnose | 尺寸预览、受控结果、可证实诊断 | 单位、测量状态、字段承诺与来源明确 |

没有“所有图层必须支持所有坐标”的承诺。不支持的组合应明确报错或标明未适用，
不能把静态源码未发现限制解释为已验证支持。不同 Position 与 link/ribbon 图元不因
名字相似强行合并；所有新增 API 仍须满足总计划的准入条件。

## 参数与返回值政策

- 保留 ggplot2 的 mapping/data/position/show.legend/inherit.aes 组合语义；
  `...` 应有明确去向、未知参数与弃用错误策略，不能藏下所有显式参数。
- 单值、按 accver 命名向量、按链设置、列表和 unit 的允许位置分别列明。不能把
  几何参数的列表语法无条件推广到 aesthetic、scale、guide 或 theme。
- 逐参数记录 NULL、NA、缺参、auto、none；不是强行规定 NULL 在全包只能有一种
  含义，而是相同职责尽量一致，合理差异须有清晰说明和错误提示。
- 公开返回结构区分文档承诺、实际返回但未承诺的附加字段、纯内部中间态。现有
  Rd 已描述的字段先按公开事实登记；收窄导出需评估迁移，不能只改分类绕过兼容。

## 内部阶段与状态的候选设计

v0.16.3 再实施显式阶段化，当前只提出以下依赖方向：

```text
原始输入/映射 -> 验证与适配 -> 全局序列/窗口布局
                         -> 每层局部几何 -> 共享障碍/标签/轨道协调
设备测量上下文 --------------------------> 受影响阶段重算
                         -> 分层输出与来源表 -> ggplot2 绘制 / 受控导出
```

抽取前逐阶段列出读取字段、写入字段、依赖、失效条件和错误出口。缓存拥有者为具体
图对象或当次 render；临时测量上下文应可嵌套还原且不残留。保留原版与候选版可比的
公开行为基线，不在同一批同时大改内部结构、数值几何与用户接口。

## 尚未确定的取舍

getter 与 export 的推荐关系仍归 v0.16.3；v0.16.2 已保留 file/files 双入口并统一解析，
保留 clean sort / prepare keep 的任务相关默认，明确来源、方向与失败行为。
v0.16.1 已决定：feature 为通用构造内核，gene 为严格定向快捷入口；保留角色美学，
固定形状使用 shape、宽度使用 width、标签使用 label_*、连接端点使用 link_*。
颜色拼写别名与重导出保留；旧 offset 和重复固定入口在 v0.16.1 删除并给迁移错误。
theme 显式 formals 聚焦专属设置，通用设置由 dots 传递并参与继承。
先看真实科研场景与已有契约，再比较保留便利入口、共享内核或渐进迁移的代价。
上述决定是本阶段可复核结果，不等同于 v1.0.0 的长期冻结承诺。

## 已补齐的共享帧依赖与状态证据

七个现有阶段的逐项输入/输出、顺序、状态归属与风险见
[audit-reviewed-details/layout-stages.csv](audit-reviewed-details/layout-stages.csv)。
这是现状合同，未实施阶段拆分：序列参考先于axis/ribbon/annotation，统一旋转先于
自动标签与裁剪，最终边界与返回装配最后执行。`compute_chord_layout()` 的形参由
被eval的表达式读取，直接函数体零引用不说明未使用；文字尺度还可能由标签阶段
更新，不应未经比较直接把各阶段改为互不影响的纯函数。

多层外部编排位于`compute_chord_geometry()`，各组保留输入、几何与子布局，之后做
track/障碍/跨层标签协调。`compute_chord_geometry_single()`内的最后一次赋值不是
“所有同类层必然覆盖”的证据；反过来，存在registry也不证明任意组合已运行隔离。

有限脚本确认context嵌套报错恢复和副本ref/scale容器隔离；源码仍有layer ggproto
引用语义及未穷举设备路径。状态审计按持久plot ref、当次geometry_cache、临时动态
render context、上游ggplot2注册/设备四类分别归属，不宣称整包没有全局副作用。

## 特征、连接和来源的现状边界

短feature可为可读性对称扩大glyph的显示角度，但输入start/end及source_row应保持；
point feature单独绘制标记。跨原点feature按circular开关拆piece，反链改变piece顺序，
不改变原始生物学链。Position先产生signed normal和lane，再由曲线局部法线或linear
坐标转换；lane号本身不带物理单位，spacing在lane_offsets阶段决定距离。

chord现状存在不同端点约定：feature/axis的角度比例用pos/lens，point link/ribbon
用(pos-1)/lens。这里只登记实际源码，未决定统一或修复；v0.16.2需对单bp、首尾和
跨原点的生物学解释提出明确方案后迁移。linear映射还有窗口source origin与x_offset，
不能只替换一个比例公式就宣称三模式合同统一。

branch在scale与after_scale后按可见样式及源端点分组；trunk保留所有source_rows，
不宣称某一个目标或统计值。当前实现依赖受支持stat/Position保持顶点行顺序；行数
改变会warning并回退，任意相同行数重排尚未得到保证。障碍provider投影实际polygons，
不将文字/leader作为gene实体；curve mapper法线正向也不恒等于相对全局原点的outward。

## 数据流程、圆图资源和输出的补充边界

read/as保留文件或对象来源，validate只报告，clean仅按显式策略处理constructor表；
prepare按固定stage顺序编排，source_rows逐阶段映回最初输入。focus/segment不拼成
连续基因组，省略间隔没有source mapping；CDS projection只映完整observed codons，
保留真实exon pieces与初始phase，不从outer envelope补造编码碱基。

圆图资源表由最终几何生成，backbone track0、physical bands、标签slots和leader
corridors属于当次build。source_row通常追输入，但某些generated片段的缺省编号仅
为build-local identity；同anchor限制位点和branch trunk另有source_rows集合。
export的original_data为所选非空组件的resolved inputs，自动装饰没有业务源行；
conversion_reports保留独立上游来源。不能把这些编号解释成统一的原文件物理行。

有限packing或leader搜索可能保留最低代价/最后候选，外部协调另写spill/collision
诊断，不能从“存在避障算法”推断全域无碰撞。元数据schema、可达性和状态生命周期
的完整来源见audit-reviewed-details；运行设备/任意staged mapping组合仍按F26等
风险定向验收，未在本版扩大用户接口或改默认图形。
