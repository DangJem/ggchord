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

本地开发仍在进行，轨道、registry 和稳定输出合同在后续批次补齐。
