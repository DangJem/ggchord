<div class="task-intro"><p>为已有输入选择局部窗口或分段视图，同时保留原始位置与显示变换的区别。</p></div>

## 本页准备 {#prepare}
在新的 R 会话中先运行下方代码，再按页面顺序继续。配色和排版是可选的显示设置，输入只选择需要的列。

{{prepare:setup,sequencepalette,docstyle}}

## 保留来源的局部聚焦 {#focus}


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
`find_ggchord_loci()` 按声明的注释/名称查询选择位置。
`link_ggchord_clusters()` 要求显式成员表和序列对，基因名称相同不构成同源证据。
`filter_ggchord_ribbons()`、`deduplicate_ggchord_ribbons()`、
`merge_ggchord_ribbons()` 整理已有比对，保留源成员和估计长度标记。
聚束和密度统计汇总已提供连接，不创建新的序列比对。

## 按源坐标放置普通注释 {#mapping}
{{code:mappingposition}}

先运行[focus 例子](views.html#focus)。辅助函数返回第一个所选窗口起点后 100 bp 位置的笛卡尔 x/y，
不把全部普通 ggplot2 坐标变成基因组 bp。分段视图应选窗口 ID，检查落在省略间隙
的位置，不跨间隙外推。区域 选择来自明确注释/同源分组 证据时见[显式同源关系](homology.html)。

## 完成后检查 {#check-result}
关注裁剪报告、断点与源坐标。被省略的间隙没有映射；显示翻转不应直接解释为生物学倒位。

<div class="next-step">继续：<a href="homology.html">连接同源基因 →</a></div>
