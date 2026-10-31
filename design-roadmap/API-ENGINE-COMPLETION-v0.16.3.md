# v0.16.3 本地开发完成记录

日期：2026-10-31。开发基线 `723575a`（v0.16.2）；实现提交为 `b8c5266`、
`337f543`、`3b5ca7c`。已完成本版全部确认范围和有限验证，按授权提交后快进整合
本地 main。未进行正式发布验收，未 push、tag、Release 或启动下一版。
GitHub 当前正式版仍为 v0.15.0。

## 实际结果

| 问题 | 完成结论 |
| --- | --- |
| F04 | 七个 quote/eval 共享帧阶段改为普通函数，显式形参和命名返回更新值上下文；保留顺序与当次参考曲线缓存；六场景同条件公开导出与原版一致 |
| F05 | 每次 layer 出现独立注册，data/mapping/固定参数和派生几何分清所有权；重复 layer、多图、复制和重复 build 已定向验证 |
| F07 | 保留 getter/export 与全部原已描述字段；推荐 export 作分层表格消费，getter 返回完整布局；无弃用或删除 |
| F18 | KEEP 既有同 ID 联合数值域/声明冲突保护、scale clone、native/mm 分离、NA/bounds/closed/baseline/分段语义；不另建轨道算法 |
| F21 | 保留圆图资源、chord外侧注释和linear共享标签的模式差异；共享低层度量/障碍/来源协议；owner、leader、track避障与显式gap已验证 |
| F22 | 四个结果入口22参数、16组返回字段明确类型、单位、来源、测量状态和缓存作用域；新增 ggchord_layout_contract 帮助主题，未增加导出函数 |
| F25 | 34项低直接调用候选联合直接/S3/ggproto/函数值/provider/字符串/quote与开发来源复核；29 KEEP、5 REMOVE；依赖语法不冒充完整运行调用图 |
| F26 | 临时动态context保留异常/嵌套还原；plot/render/ref/scale/layer归属明确；所有布局检查设备/context失效；预览/诊断副本保护源缓存 |

F48 修复 `build(p)` 后 export 序列为空、复用 layer 共 ID 和派生状态写回；基线教学图
分别为先 build 后 export 0 行、直接 export 500 行，当前两路径均保留序列几何。
F49 修复 preview 覆盖原 plot 缓存；F50 为缺失标签的空框加保护，消除 Inf/-Inf 警告，
可见几何和候选顺序不变。F48 仅增加一个必要的公开 testthat 回归；本轮未执行 testthat。

五个删除项均为无调用私有 helper：`key_glyph_gene`、`ggchord_enforce_label_side`、
`ggchord_is_identity_position`、`ggchord_restriction_label_lanes_cartesian`、
`ggchord_preview_height`。保留9个S3方法、.onLoad、图例/阶段/provider回调、REBASE
数据生成器、label_box_conflicts检测器和有实际调用的if_null_else。无公开删除版本。

v0.16.1 的固定样式优先、同类层独立、手动 scale 优先和 staged mapping 保持；
v0.16.2 的原 ID/源行、未知链、单bp方向和聚合估计标记未改。历史 chord 锚点差异、
chord NULL原点控制与已正确默认外观保持，没有将内部整理改写成生物学坐标迁移。

## 实际验证与重现

环境：R 4.6.0（2026-04-24）、ggplot2 4.0.3、roxygen2 8.0.0。

| 验证 | 结果与实际边界 |
| --- | --- |
| `Rscript tools/validate-engine.R capture <临时.rds>`，随后 `compare <临时.rds>` | 六场景PASS；capture从723575a的原阶段直接重建，无第二checkout；两批阶段抽取及最终比较保留公开导出。相同seed=163、8×6英寸、Helvetica和PDF度量设备；覆盖默认chord、正负curvature/非单位radius/反向orientation、radial/auto、圆图资源、linear和分段track。临时浮点快照不进仓库 |
| `Rscript tools/validate-engine-state.R` | 八项PASS；构建后导出/复用层、sibling/render copy、restriction派生参数、嵌套错误context、两个PNG目标尺寸（6×4和9×5）物理4mm band、显式静态preview、诊断失败设备恢复、普通标签cache-only失效。不是发布设备矩阵 |
| `Rscript tools/validate-tracks.R` | 既有三坐标/分段/NA与五项新风险PASS；同ID共享domain/source layers、声明冲突、log10/OOB/baseline、clip/NA/原坐标、完整/缺测closed seam。没有穷举参数等价格式 |
| `Rscript tools/validate-engine-registry.R` | 四项PASS；两种坐标独立owner/source_row、linear placed/hidden与leader证据、track band自动避障/显式gap优先、post-scale trunk来源集合与变化staged样式 |
| `Rscript tools/validate-fixed-styles.R` | PASS；重建layer后显式固定值/默认映射、feature wrapper、继承与after_scale保持 |
| roxygen/Rd与源码 | 仅roxygen Rd/namespace生成；85份R与92份Rd无警告解析，diff无空白错误；未重建冻结网站 |
| `Rscript tools/audit-engine.R` | 七阶段、187项输入字段来源、1,224条具名内部依赖语法、34项候选、22参数、16组返回字段与来源MD5；无未决删除候选，静态生成器不生成运行PASS |

基线脚本的四个反向ribbon输入验证警告为预期，保留其真实方向；未把警告当作构建失败，
也未为得到无警告输出改原数据。缺失标签的空框警告修复后不再出现。

机读证据为 `audit-generated/engine-*-v0.16.3.csv`，人工字段/去留依据为
`audit-reviewed-details/engine-*-v0.16.3.csv`。最终结构快照基于实现提交`3b5ca7c`，
逐文件散列识别实际源码；旧`api_*.csv`、grammar/data表保留各自历史来源，不冒充新版证据。
详细契约与生命周期见 [引擎与输出契约](API-ENGINE-CONTRACTS-v0.16.3.md)。

## 已明确的边界与下一版建议

本版无未完成开发项或待用户选择的兼容冲突。现有公共返回字段全部保留；新增附加列
不意味着固定列数合同，resource ID/source_row不是永久生物ID/统一文件行号。
测量panel与geometry收敛分开；普通R赋值保留ggproto引用语义，用户任意手动修改
共享layer环境或quosure产生副作用不在自动隔离保证内。

有限场景无已知严重状态串扰，不代表任意几何、字体、平台或设备穷举。未观察到需要
定向算法优化的实测退化，因此未作无证据的性能改写、固定秒数承诺或极端压力测试。
候选阶段仍需按实际风险处理复杂布局/平台兼容；正式发布须用户明确授权后运行完整
testthat、R CMD check、英语/汉语vignette及PNG/PDF/SVG验收。本轮均未运行。

建议**无需独立启用v0.16.4**：前版参数迁移已有执行方向，本版getter/export字段和
缓存边界已收束，没有足以支撑独立稳定化版本的已确认剩余开发；v0.16.5承担既定的
教程/README/网站/图片源重构，v0.16.6承担候选审计和完整发布前收敛。
这是建议，是否跳过仍由用户决定；保留v0.16.4空号及后续编号，不自动启动任何下一版。
