# v0.16.3 引擎与输出契约

开发基线：`723575a`（v0.16.2 本地完成）。本文件区分源码合同、有限运行证据和
后续正式发布验收；旧 `api_*.csv` 仍是 v0.16.0 快照。

## 显式阶段与依赖

七阶段均为普通函数，形参是读取合同，返回的命名 list 是唯一的跨阶段更新。
`compute_chord_layout()` 将已规范化参数放入当次调用拥有的值上下文，按函数形参
选择参数，以返回值更新上下文；不向函数帧注入变量，不执行 quote 阶段。
循环变量、采样候选和局部闭包留在所属阶段，不能偶然覆写后续状态。

```text
规范参数 ─> sequence ─> axis
                 ├───> ribbon
                 └───> annotation
sequence/axis/ribbon/annotation ─> transform ─> label ─> finalize
                                            ↑          ↑
                                设备文字度量/固定障碍   最终文字尺度与可见几何
```

逐字段读取/写入见 `audit-generated/engine-stages-v0.16.3.csv`。
顺序保留 sequence→axis→ribbon→annotation→transform→label→finalize；
axis/ribbon/annotation 共享 sequence 结果，但没有相互几何依赖，不据顺序捏造依赖。
transform 按紧凑几何估计文字尺度，label 可更新尺度，finalize 使用更新后的值。

| 阶段 | 返回与所有权 | 单位/失效条件 |
| --- | --- | --- |
| sequence | 起止角、参考曲线、弧、局部放置闭包、转换后的轴距离 | bp 长度、弧度、布局坐标；参数变化重建；参考缓存只属于外层当次 geometry 调用 |
| axis | 轴线、刻度与文字锚 | 布局坐标、角度；序列参考、break、样式或物理测量变化重算 |
| ribbon | 多边形与描边/线型使用标记 | 布局坐标；端点、方向、障碍和控制点改变重算 |
| annotation | region/feature 多边形、标签锚、色板与序列文字 | 原 bp 和来源列保持；Position 与物理文字度量影响显示 |
| transform | 旋转后几何、紧凑范围和文字尺度 | 布局坐标、度、布局单位/英寸；不改原 bp |
| label | 可见标签、leader、裁剪尺度与 fit 诊断 | 布局坐标和 mm 字号；尺寸、字体、障碍或策略改变重算 |
| finalize | chord_layout | 汇总边界与既有返回成员；不把候选搜索状态导为稳定承诺 |

## 可比基线

`tools/validate-engine.R capture <临时.rds>` 在抽取前保存六场景公开 export；
`compare` 在相同 seed=163、8×6 英寸、Helvetica 与同一 PDF 度量设备上比较。
临时 RDS 不进入仓库，不把浮点或文字框结果固化成永久测试。
第一批 sequence/axis/ribbon 与第二批其余阶段后均通过六场景一致比较。
输入为脚本内明确教学数据，包括未知链、不同映射、反向端点、同类独立层、
不同 radius/curvature/rotation、圆图资源与分段 track；不是所有几何/设备穷举。

## 已复现并修复的状态问题

- F48：基础图先 build 再 export 返回 0 行序列，直接 export 返回 500 行。
  layer ID 只在部分构建副本赋值，出口又依赖原图层 ID，导致调用顺序改变结果。
- 同一路径还会给用户 layer ggproto 写入 ID 和 restriction 的派生参数。
  已改为每次构建局部重建 layer，并用一致的注册顺序供 build/export 使用。
- F49：明确尺寸 preview 会覆盖原图缓存，返回原设备后 cache-only getter 失败；
  已改为渲染副本。所有布局在设备/context 变化后失效，普通标签不再跳过检查。

`tools/validate-engine-state.R` 的八项有限检查通过：构建后导出/复用 layer、
多图/副本、restriction 参数、context 嵌套报错、两个 PNG 尺寸的物理轨道、
静态预览、设备构建失败和普通标签缓存失效。原图的业务表、mapping、固定参数、
图层参数与缓存由对应检查保留；不是任意用户 quosure 都无副作用的保证。

## 状态所有权与生命周期

| 对象 | 拥有者与有效期 | 更新规则 |
| --- | --- | --- |
| 用户 data/mapping/参数 | 每个原始 layer；constructor 的三表仅作对应 fallback | 构建时解析副本；同一 layer 重复加入也独立注册，不写原 layer |
| layer ID/几何 registry/子布局 | 当前 plot 规格的一次 geometry 计算 | 确定性按图层位置分配；export 在局部副本使用同一规则；ID 不是跨编辑的永久主键 |
| plot ref | constructor/新组合拥有；普通 R 赋值可共享引用 | 普通 build/getter 可填充缓存；新组合和 render copy 分配新 ref |
| scale 训练 | 几何计算或 ggplot build 副本 | 手动 scale 优先；位置/shape/track scale 显式 clone 后训练；原 scale 不被训练 |
| reference cache | 一次外层 compute_chord_geometry | 多个子布局可复用同一曲线 key；调用结束不跨图保留 |
| render context | 动态调用栈临时持有 | on.exit 嵌套/异常还原；不是全局布局缓存 |
| 设备/gtable 度量 | 一次目标 render | 尺寸/context 改变使估计失效；diagnose/preview 使用独立副本，不覆盖源缓存 |

普通 R 赋值保留 ggplot2 的 ggproto 引用语义，不能承诺用户直接修改共享
`p$layers[[i]]` 环境后，另一个赋值副本完全不受影响。支持的组合 `q <- p + ...`、
重复 build、内部 render copy 的隔离已经定向验证；不宣称任意手工内部突变有自动追踪。

## 轨道、registry 与共享边界

F18 保留既有轨道实现：同 track ID 的 side/order/width/gap 必须一致，数值域跨
层/序列/窗口联合训练并包含 bar baseline；数值 scale clone 后训练，limits 记录原值域。
物理 mm 与 native 布局距离分开保存，value 域不决定 band 的空间位置。
`tools/validate-tracks.R` 复用既有三坐标/分段/NA 场景并补五项风险检查：同 ID 域、
声明冲突、log10/OOB/baseline、bounds/NA/来源、closed/NA seam。全部通过。
没有据“已有保护”推断任意组合无缺陷，也没有另建一套轨道算法。

F21 保留三种模式的显式差异：圆图 registry 分配内侧物理带与外侧 slot/corridor；
chord 共享外侧注释障碍；linear 在水平骨架/窗口 frame 内协调固定/自动文字与 leader。
共享的是文本度量、框/线相交和原始来源协议，不强行合并坐标专属资源或候选算法。
`tools/validate-engine-registry.R` 核查独立层重复 source_row 的 owner、圆图资源、
线性 placed/hidden/leader 诊断、track band 自动避障与显式 gap 优先，以及分支的
源行集合/可见 staged 样式；四项通过。
缺失文字会产生空候选框，旧 min/max 触发 Inf/-Inf 警告；现加空框保护，六场景
公开导出与基线一致。没有改变可见标签的候选顺序或布局数值。

## getter、export、preview 与 diagnose 的结果合同

保留四个入口和当前签名；22 个参数的当前 formals/旧语义来源与所有权另存
`engine-parameters-v0.16.3.csv`，16 组返回字段的类型/单位/来源/作用域见
`engine-return-contracts-v0.16.3.csv`。完整用户帮助为 `ggchord_layout_contract`。
推荐表格消费用 export，完整布局用 getter；不是弃用/收窄，没有公开删除版本。
已有 Rd 中 annotation_registry、sequence_layout、轨道、诊断和来源字段全部保留。

- 几何与 bbox 是布局 data 单位，原 bp 不重编码；文字与物理贡献为 mm，角度为度。
  metadata 明确已应用 layout rotation、尚未应用 ggplot2 coord transform。
- source_row 是该层 resolved input 的位置；源文件物理行、GenBank record/feature、
  聚合 source_rows 与自动装饰的 NA 来源分别解释，不伪称统一原文件行号。
- original_data 只返回所选非空组件的 resolved inputs，含角色映射/适配；conversion
  属性保留上游对象 provenance，不是原始业务文件归档。
- 物理 track registry 的 `_mm`/native/normal 字段和数值 limits/unit 各有意义；
  圆图 track/slot/sector 是 build-local 资源号，不是命名定量 track ID 或生物 ID。
- render context 的 panel status 与 geometry_metric_status/converged 分开；实测 panel
  不意味着字体嵌入、任意 grob、跨软件透明度或所有碰撞均已通过。迭代次数/私有缓存
  不作为稳定返回合同；缺测字段与空表保持既有语义，不用固定列数取数据。
- preview 是临时静态图片路径；diagnose 是单个目标设备的有限证据；正式保存仍用 ggsave。

## F25 内部候选与性能政策

当前源码 AST 产生 1,224 条去重的具名内部依赖语法和 34 个低直接调用候选。
联合 NAMESPACE S3、ggproto/provider/函数值回调、字符串与引用表达式，以及
`tools/`、`tests/`、`data-raw/`、示例/帮助证据复核，29 项 KEEP、5 项 REMOVE。
完整人工理由在 `audit-reviewed-details/engine-unused-review-v0.16.3.csv`，生成表
在 `engine-unused-candidates-v0.16.3.csv`；依赖图仍标 SYNTAX_ONLY，不冒充完整运行图。

删除的是 key_glyph_gene、ggchord_enforce_label_side、ggchord_is_identity_position、
ggchord_restriction_label_lanes_cartesian、ggchord_preview_height 五个私有无调用 helper。
保留 .onLoad、9 个 S3 方法、stage/legend/provider 回调、REBASE 数据生成器、验证器的
label_box_conflicts 和确有直接调用的 if_null_else；不是按零词法计数自动删除。
没有新增公开函数/参数或删除原有导出；帮助主题不增加导出名称。

未观察到需要几何算法优化的实测性能/内存退化，本版不声称加速或省内存，也没有
极端压力/固定秒数目标。复杂几何质量、平台、字体和完整设备矩阵继续留给显式
候选/正式发布验收；不把本次有限正常输入当作性能基准或全域几何保证。
