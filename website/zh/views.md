## 保留来源的局部聚焦 {#focus}
先运行[共享准备代码](index.html#quick-start)：

{{code:focus}}
{{figure:focus}}

focus 把保留位置移到局部一基坐标，同时保留源原点、截断标记与报告。
`boundary = "trim"` 可按线性显示投影裁剪连接，不重新比对序列；缩放后的比对长度
带估计标记。`coordinates = "source"` 恢复源坐标轴标签，重复聚焦后同样适用。

## 独立窗口与局部翻转 {#segments}
{{code:segment}}

各窗口有可逆的源映射，被省略的间隙没有映射。`map_ggchord_position()` 把显式源
位置映射到图中，供普通 ggplot2 注释使用。多个窗口共享源序列时应指定窗口 ID。
显示方向不表示原始链发生改变。

## 显式生物学关系 {#homology}
`find_ggchord_loci()` 按声明的 feature/名称查询选择位置。
`link_ggchord_clusters()` 要求显式成员表和序列对，基因名称相同不构成同源证据。
`filter_ggchord_ribbons()`、`deduplicate_ggchord_ribbons()`、
`merge_ggchord_ribbons()` 整理已有比对，保留源成员和估计长度标记。
聚束和密度统计汇总已提供连接，不创建新的序列比对。
