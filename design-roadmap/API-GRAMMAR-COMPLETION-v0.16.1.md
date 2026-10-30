# v0.16.1 本地开发完成记录

日期：2026-10-31。状态：已完成已确认的开发范围和相关轻量验证；未执行正式发布
验收，未发布。GitHub 当前正式版仍为 v0.15.0。按用户指示，完成提交整合本地 main，
不自动推送、发布或启动下一版。

## 实际结果

- feature 是通用构造内核；gene 是其定向箭头快捷入口，保留严格正负链、gene_fill
  独立美学和原默认几何。plasmid 保留有独立默认语义的便利组合。
- 固定形状唯一使用 shape，gene/feature 宽度统一 width，gene/feature 标签公开
  参数统一 label_*，line/ribbon 的端点参数统一 link_gap/link_ctrl_point。
- 在 v0.16.1 删除旧固定 feature_shape、feature_width/gene_width、gene_label_*
  几何参数、ribbon_gap/ribbon_ctrl_point、gene_offset/feature_offset 及 trunk 兼容值。
  旧调用明确给迁移错误；删除旧 offset 的内部解析路径。用户已取消旧 v1.0.0 期限。
- 显式固定样式优先，省略默认时保留映射；正确处理局部/继承映射、wrapper、
  scale 推断和 after_scale。固定序列 colour/color 作用于 seq_colour。
- restriction straight/elbow 真正分派为直线/单折点，内侧 elbow 向内离开；radial
  默认与旧源码对照一致。直线不自动避障，不承诺全域无碰撞。
- theme 自身签名只列专属设置，通用 ggplot2 命名设置通过 dots 继续生效并参与
  角色继承。补全插件可能自行展开上游 dots，未认证各 IDE 插件行为。
- 几何角色拒绝 late stage 和命名空间形式的 after_stat/after_scale；视觉 staged
  mapping 仍按 ggplot2 求值。澄清单位、职责、序列轴与独立 bp 参照条、重导出来源。
- 保留英美拼写别名与22个有独立工具职责的重导出，不按数量裁剪入口。

## 证据与重现入口

环境：R 4.6.0（2026-04-24）、ggplot2 4.0.3、roxygen2 8.0.0。

| 验证 | 实际结果与边界 |
| --- | --- |
| `Rscript tools/validate-grammar.R` | 14项 PASS；三坐标独立双 gene 层与顺序、feature/固定/repel 标签组合、chord/collinear staged density、F39 分派与内侧反向/开口、F40 标量域。教学重叠连接产生2个预期验证警告，不是构建失败 |
| `Rscript tools/validate-fixed-styles.R` | PASS；固定与省略值、继承 link 映射、feature/plasmid 显式形状、seq/gene/region 样式、after_scale 优先级；无默认 scale 警告 |
| 变更前源码直接对照 | 从 `9c2b5fa` 提取 theme/gene 定义，5个主题默认和通用文字/图例覆盖一致，gene 默认布局一致；radial 路由与变更前源码一致。只比较本次代表输入，不把具体坐标固化为测试 |
| 迁移与错误直接脚本 | 新版完整迁移示例 build 通过；旧 offset/固定形状/标签/连接参数得到替代提示；gene 未知链仍拒绝；命名空间 late stage 角色明确拒绝 |
| Rd/R 源码解析 | 变更的 Rd 和全部 R 源码解析通过；只运行 roxygen Rd/namespace 生成，无包构建或发布检测 |
| `Rscript tools/audit-grammar.R` | 72个绘图/构造导出、1,038行当前参数职责、32项美学/角色、69项三坐标组合；输入与生成器散列核对一致 |

本版产物位于 `audit-generated/grammar-*-v0.16.1.csv`。兼容矩阵逐项区分源码合同与
有限运行组合，未运行项明确标记 NOT_RUN。旧 `api_*.csv` 是 v0.16.0 审计快照，
不把旧版签名或行号冒充重构后的当前证据。

最小回归复用现有 ribbon 行为案例补充显式 alpha 优先，并增加一个 leader 选择
行为案例；其余测试只迁移调用。没有运行 testthat、devtools::test、R CMD build/check、
vignette 重建或 PNG/PDF/SVG 发布设备矩阵，也未重建冻结 docs 网站。

## 逐项完成与后续边界

F06、F08–F15、F19、F27–F29、F39/F40 均有明确实施或 KEEP 结论；矩阵、职责、
NEWS、roxygen/Rd、迁移、版本及决策记录已同步，无尚待用户取舍的本版阻塞项。

本地开发完成不等于发布验收通过或已发布。数据读取/生物学契约仍归 v0.16.2，
内部生命周期、缓存、任意 stat/Position 组合及输出字段收敛归 v0.16.3，网站/图片源
重建归 v0.16.5，平台/完整设备收敛归 v0.16.6。现存历史 README 图片生成器使用
更早的已退役 API，已明确标记为历史资产溯源代码，不冒充当前可运行教程。
