## 职责分配 {#grammar}
| 任务 | 接口 |
| --- | --- |
| 把数据列映射到 aesthetic | `aes()` |
| 定义映射规则与图例 | `scale_*()`, `guides()` |
| 几何与固定样式 | `geom_*()` |
| 局部放置 | `position_*()` |
| 全局坐标系 | `coord_*()` |
| 非数据外观 | `theme_*()`, `theme()` |
| 标题与说明 | `labs()` |

固定样式可直接传给 geom，数据驱动差异则使用映射和 scale。显式固定值优先于映射，
省略固定参数则保留映射；用户 scale 优先于自动默认值。
`shape` 是固定 feature 几何参数，`feature_shape` 仍是映射角色。
显式 `shape = NULL` 无效，不表示改用映射。

## 映射连接相似度 {#mapping}
先运行[共享准备代码](index.html#quick-start)：

{{code:adjustments}}

## 放置 feature {#position}
使用 `position_strand(c("+" = .1, "-" = -.2))` 指定带符号的逐链位移。
`position_feature_stack()` 处理重叠 feature，`position_track()` 分配命名定量轨道。
逐序列几何参数可用单值、等长无名向量或按 `accver` 命名的向量；推荐命名向量。
这套语法不扩展到 scale、guide 或 theme。

每层独立拥有 `data`、`mapping` 和参数，同类 feature/连接层不会互相替换。
`show.legend` 和 `inherit.aes` 遵循 ggplot2。输入角色映射在几何生成前处理，
后期视觉映射遵循 ggplot2 求值。详见 [API 契约](reference/ggchord_api_contract.html)。
