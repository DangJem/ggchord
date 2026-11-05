<div class="task-intro"><p>减少连接拥挤，同时保留每一步的来源：筛选、去重、合并和聚束有不同含义。</p></div>

<span id="prepare"></span>

## 从原始连接表开始 {#input}

{{code:ribboninputs,ribbons}}

{{results:ribbons}}

前三列沿同一条处理路径统计。聚束从原始表开始，行数单独理解。

## 理解每步转换 {#transformations}
筛选保留满足已声明阈值的观测；互惠重叠去重选择已有代表观测，不累加相似度；
合并把邻近兼容区间组合为显示外包范围；聚束按邻近 query/subject 位置形成汇总。
逐步检查 `$report` 和源行成员，特别注意冲突的额外元数据。
合并或聚束范围不是重新计算的核苷酸比对，估计长度保留标记，不凭空生成相似度。

## 密度是显式统计 {#density}
{{code:ribbondensity}}

{{figure:density}}

这里把两端各分成四个位置箱。共享同一序列对、方向和位置箱的观测越多，带越深。图例显示每箱观测数；`after_stat(bundle_n)` 取统计层算出的计数。几何仍是一行输入对应一条带。需要把邻近带合成一条时，用 `stat_ribbon_bundle()`。

## 选择区间区间连接或点连接 {#points}

{{code:pointlinks}}

{{figure:pointlinks}}

点连接只画前三条真实比对的中点，箭头指向目标端。需要显示整个比对区间时，继续用连接带。

## 完成后检查 {#check-result}
先检查各步骤 report 和成员行，再读图。聚合带的范围不是新比对；点连接是已有 HSP 中点的摘要。

<div class="next-step">继续：<a href="views.html">查看局部和分段区域 →</a></div>
