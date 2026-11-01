# v0.16.6 当前合同审计与取舍

基线为正式版 v0.16.5 重发标签 `d697345`，候选比较基线为 v0.16.4。当前开发包
为 `0.16.6.9000`，并非已发布 v0.16.6。本文的审计结论不授予发布权限。

## 覆盖与证据等级

`tools/audit-rc.R` 从当前 namespace、formals、Rd、NAMESPACE 和源码位置生成
`audit-generated/rc-*-v0.16.6.csv`，不改历史候选。130 个导出（108 自有函数、
22 上游重导出）、1,236 个自有形参、49 项美学角色、9 个 S3、130 返回说明及
200 迁移记录均有关联合同。形参数量只是覆盖度，不是参数全排列执行证明。

- `rc-exports/parameters/s3`：当前签名、默认、所有权、别名和文档的静态/对象身份
  核对。`rc-parameters` 保留原语义表键及当前 Rd 文案，逐项呈现差异。
- `rc-aesthetics`：保留输入、视觉及计算后角色边界。源码字面量存在是定位检查，
  不是每个入口都支持该 aesthetic 的证明。`R/zzz.R` 的 mapping 规范化与各 geom
  的角色列表、`R/stat-ribbon.R` 的 computed 字段、`R/geom-track.R` 的 scoped
  track key 分别复核；私有 key、annotation 元数据和 draw 坐标不新增为公开 aesthetic。
- `rc-returns` 记录当前 Rd；`rc-data-fields` 的 34 个入口合同和
  `rc-engine-fields` 的 16 组字段继续关联 v0.16.2/3 的单位、来源、空值及缓存合同。
  公开 annotation_registry、sequence_layout、conversion、track 与诊断字段保留。
- `rc-migrations` 逐条继承当前决定，保留首删版本与“最迟已删除”上界的区别。
  本版没有新增重命名、弃用、删除或同义入口；只有代表 guard 实际触发。
- `rc-compatibility` 的旧支持/运行证据列完全保留；新 `rc_evidence` 仅关联确实
  出现在六组运行中的 geom/coord 对，其他项仍为 NOT_RUN_IN_RC。
- `rc-behavior`、`rc-format-validation` 是当前有限运行；`rc-support` 区分声明、
  安装情况、实测与未知。源码 MD5 是字节指纹，不证明语义正确。

## 与候选之间的实际差异

签名、5 个直接拼写别名、22 个原样上游绑定与 9 个 S3 注册均无漂移。v0.16.5
相对候选的参数/返回文案变化主要为输入最小列、标准 BLAST、来源字段、主题物理
单位、arrow/Position、识别与投影能力限制、固定/映射优先和上游委派说明；
`rc-parameters/returns` 的 changed 列逐项保留，不将文案差异当作新实现。

必须另行追踪的 v0.16.5 行为事实：

| 项目 | 当前源码与合同结论 | 证据 |
| --- | --- | --- |
| F52 零行引物组合 | input transform 的零行保护；build-local 空层不求值缺失映射，原层与导出登记保留 | R/geom-primer.R、R/ggchord-frame.R；既有单一回归；v0.16.5 正式验收为历史运行 |
| F53 四列 gene | 缺失 anno 只清 annotation-based feature 默认映射，strand-based gene_fill 保留 | R/zzz.R；本版 layers 实际验证手工 +/- 标度 |
| 稀疏混合引线 | 有外侧酶切与外侧 primer 时即共用真实锚点排序；不再要求 30 个位点 | R/layout-external-annotation.R；本版 circular 完整 feature+primer+site 无引线交叉 |
| export 声明 | geom_feature_label_repel 的 roxygen @export 恢复，当前 NAMESPACE 及帮助完整 | R/geom-feature-label.R、NAMESPACE；当前导出身份/文档检查 |
| 帮助编码 | ASCII 显示副本替换 α 仅用于默认 PDF 帮助例子，原记录不改 | R/api-contract.R 等当前示例；v0.16.5 远端事实见正式验收 |
| utils 限定 | ribbon smooth refinement 中 head/tail 明确使用 utils | R/layout-ribbon-step.R；没有新 API/默认变化 |

旧发现表中的“待选择”按 v0.16.1–4 后续决定收束，不把旧建议恢复成当前任务。
F52/F53 不重新编号或冒称本版首次修复。

## 保留入口的理由

| 当前入口/角色 | 取舍与职责 | 源码/原决定 |
| --- | --- | --- |
| get_chord_layout / export_ggchord_layout | KEEP；完整布局与按图层表格消费不同。表格推荐 export，不弃用 getter 或删除其成员 | F07；R/zzz.R、R/layout-export.R、ggchord_layout_contract |
| view_ggchord / diagnose_ggchord | KEEP；临时静态图路径与目标尺寸有限测量证据不同；正式保存仍用 ggsave | F22/F26；R/viewer.R、R/diagnostics.R、R/render-context.R |
| gene / feature / plasmid 便利入口 | KEEP；gene 严格定向，feature 保留未知链，plasmid 具有已说明的 width/Position 默认。共享几何不抹平输入及预设差异 | F10–F13；R/geom-gene.R、R/geom-feature.R |
| point link / ribbon | KEEP；点连接与成对区间面积不同，avoid 默认分别 none/uniform，NULL chord 控制点仍为原点 | R/geom-link-line.R、R/geom-ribbon.R、API-CANDIDATE-v0.16.4.md |
| 5 colour/color 直接别名 | KEEP；原样函数身份，沿用 ggplot2 拼写习惯。文档优先 colour；不是新 API 冗余迁移 | F14；seq manual、ribbon manual、link manual/gradientn、guide colourbar；rc-exports 全部 identity 核对 |
| 22 重导出 | KEEP；20 ggplot2 + 2 grid 的常用组合工具，用户无需重复 attach。原始上游对象 identity 全部核对；不冻结上游未来签名 | F15；R/reexports.R，rc-exports ownership 逐项列出 |
| S3/ggproto | 9 注册 KEEP；泛型派发为公开组合协议。geom、coord、Position、Stat/provider 回调内部字段不承诺永久结构 | NAMESPACE、rc-s3、R/aaa-geom.R、R/stat-ribbon.R、R/coord-*.R |

## 灵活格式的逐参数边界

27 条当前公开参数组合见 `rc-parameter-formats-v0.16.6.csv`。保留所有既有有效
格式，没有删除数量目标，没有已决定的弃用/移除日期。推荐按 accver 命名向量：

- sequence helper 支持单值、按当前显示顺序的等长无名向量、ID 命名向量、
  ID/序号命名列表及无名列表。ID 匹配优先于序号；`list(value)` 循环。
- gene width/固定 label 四参数增加 +/- 链向量；逐序列列表元素为链值对，链值对
  内部可用无名二元向量。不同于 `list(value)`，`list(A=value)` 不是这一 helper
  的逐序列链值对格式。未指定键保留参数 fallback，不能统一解释 NULL 为关闭。
- position_strand 的单值按 +/- 展开；显式链值本身带符号，缺省链/序列补零。
  逐序列列表还允许单值；嵌套链值对必须具名。position_plasmid 是共享位移，
  拒绝链值对。不能把 width helper 的全部等价语法推广给 Position。
- position_track 的数值 width/gap 使用 sequence helper，但缺少任何显示序列值
  会拒绝；物理 unit 只能是一个绝对尺寸。其数值域不决定轨道物理位置。
- 轴和 theme、aesthetic、scale、guide 不接收上述扩展。轴 theme 的 orientation
  当前只允许单值；旧私有 process_axis_orientation 的向量分支不是公开合同。
- 名称应唯一。纯数字 ID 与序号列表、ID 为 +/- 与顶层链简写可能有歧义；优先
  使用明确 ID 命名的嵌套结构。重复键不是已承诺的独立多值能力；本版不静默
  统一历史 helper 的覆盖次序，也不新造重复键语义。

`validate-rc-formats.R` 只运行一个真实基因数据、命名半径重排加嵌套链值的公开组合，
不按每个 helper 或格式穷举，也不把静态表改为逐格式 PASS。

## 职责与所有权复核

aes 输入角色先于几何求值，视觉映射可在已支持阶段使用 after_stat/after_scale；
Stat 只消费已计算 ribbon 统计字段，未持有设备碰撞状态。Position 声明序列局部
位移/逻辑 lane/轨道，coord 管物理空间及共享障碍，scale 管值域/映射/guide，
theme 管固定非数据外观。注册表保留 documented 诊断结果，不公开候选求解状态。

`compute_chord_geometry` 的每层注册与 build-local clone 继续保证独立 data/mapping/
参数；用户 scale 在 clone 中训练且优先于自动默认。相同 layer 对象重复加入也各自
注册。支持的组合与 render copy 隔离，不承诺手工突变共享 ggproto 环境自动隔离。
本版 layers/output 有实际来源集合、用户颜色、原对象/缓存及报错设备还原的证据。

来源 source_row 是 resolved input 的行位置；原文件行、source_rows 聚合、NA 装饰
另有意义。生物 strand、显示 orientation、link direction 不合并。CDS 邻域与显示
投影不是重新计算比对；残基投影不是原生 aa 绘图。

## 本版问题与风险状态

| 编号 | 类型 | 结论 |
| --- | --- | --- |
| F54 | 确证 P2 正确性缺陷 | 已修复 inside primer 被外侧分配器移动；一个公开回归，六组内 circular 有实际断言 |
| F55 | 确证 P2 诊断缺陷 | 已修复 inside restriction registry 的 outer 误标；复用酶切案例补侧别断言，circular 实际验证 |
| RC-R01 | 参数匹配风险 | 已明确位置写法随当前显示顺序，特殊键及重复键边界；推荐命名向量，不新承诺重复键行为 |
| RC-R02 | 支持证据限制 | 下限版本、其他 OS/设备/字体、可选依赖路径未知，单独 NOT_RUN；不以此认证本版发布，也不伪造成功 |
| RC-R04 | 已观察的非阻断视觉限制 | 完整 pUC19 的显式 inside 引物长文本与部分 feature 区域交叠；侧别/内容正确，默认 max_overlaps=Inf 不保证跨 feature 的零交叠。推荐默认 outside 或减少显示文字。当前有限场景不认证 inside 全局避障；进一步全局内侧避障属于待确认优化，不承诺下版 |
| RC-T01 | 本轮审计工具问题 | 已修正 validate-rc.R 的报告闭包读取初始结果而写 NOT_RUN；最终文件必须由显式传入实际 result 的运行生成，没有手工批量升级 PASS |
| RC-R03 | 有限几何限制 | 已检查稀疏但局部拥挤的完整 pUC19、三标签和同类多层；密集最多40真实cuts为显式工具场景，尚无全域无交叉保证 |

没有未解决的已确证严重正确性问题或需要用户选择的破坏性方案。未知支持组合是
明确的发布验证边界，不表示已验收；新发现的确定问题仍须独立处理。

## 有限视觉/性能工具

旧脚本的写入、副作用与选用结论见 `rc-tools-review-v0.16.6.csv`，没有批量执行 tools。
新工具始终写 rc 版本文件或忽略目录，历史 CSV 与 docs 不改。

- `Rscript tools/render-rc.R circular`：单张 PNG，9×8 in、120 dpi、sans；已检查
  完整文本、标签/引线、方向和图边界。inside 也已渲染，侧别/内容通过，跨 feature
  交叠按 RC-R04 如实记录；dense/chord/collinear/track 可显式选用；
  每次只选一个场景。dense 为真实非冗余单切酶显示子集，最多40行，不绘制全部命中。
- `Rscript tools/benchmark-rc.R 1`：真实13个 HSP、完整圆图、真实 AmpR 邻域及100 bp
  GC bins 各一次 build；8×8 in Helvetica PDF 只作度量。输出时间与构建对象占用，
  后者不是峰值内存；没有跨机器固定秒数门槛、加速声明或旧版本性能比较。
- `Rscript tools/validate-rc.R`：六组有限行为；不等于 testthat 套件、vignette 或设备矩阵。

完整发布验证仍需另行授权；本地完成及提交依据见 `API-RC-COMPLETION-v0.16.6.md`。
