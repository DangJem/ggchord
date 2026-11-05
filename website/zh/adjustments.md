<div class="task-intro"><p>把常见调整放到正确的 ggplot2 组件：数据映射、几何位置、坐标和主题各负其责。</p></div>

<span id="prepare"></span>

## 职责分配 {#grammar}
[绘图基本概念](modes.html#prepare)解释图层、映射、标度和坐标的分工。下面只演示如何把这种分工用到比对相似度。

## 映射连接相似度 {#mapping}


{{code:adjustments}}

{{figure:mappingexample}}

## 放置注释 {#position}
使用 `position_strand(c("+" = .1, "-" = -.2))` 指定带符号的逐链位移。
`position_feature_stack()` 处理重叠注释，`position_track()` 分配命名定量轨道。
逐序列几何参数可用单值、等长无名向量或按 `accver` 命名的向量；推荐命名向量。
这套语法不扩展到标度、图例或 theme。

每层独立拥有 `data`、`mapping` 和参数，同类注释/连接层不会互相替换。
`show.legend` 和 `inherit.aes` 遵循 ggplot2。输入角色映射在几何生成前处理，
后期视觉映射遵循 ggplot2 求值。详见 [接口说明](reference/ggchord_api_contract.html)。

## 完成后检查 {#check-result}
旋转整体图形使用 coord；逐序列几何使用明确命名值；颜色随数据变化时使用映射和标度。

<div class="next-step">继续：<a href="parameters.html">灵活参数配置 →</a></div>
