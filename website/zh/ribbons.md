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

点与区间连接的几何选择见[连接课](multi-sequence-links.html#interval-points)。

## 完成后检查 {#check-result}
先检查各步骤 report 和成员行，再读图。聚合带的范围不是新比对；点连接是已有 HSP 中点的摘要。

<div class="next-step">继续：<a href="views.html">查看局部和分段区域 →</a></div>

需要明确理解 bundle_n、bundle_weight、density 与 bins/weight/group_by/min_bundle，继续[连接统计](ribbons.html#statistics)。

密度着色与聚束不是同一个操作。先保留原始 HSP，再说明箱数、权重、分组与显示汇总的意义。

## 真实输入与最小字段 {#statistics}

本页独立使用随包的完整序列长度和原始 HSP，不运行新比对。必需 ID 与两端坐标之外，weight=length 需要 length，weight=pident 同时需要 length 与 pident，不能从端点猜这些统计量。

{{code:comparedensitydata}}


## 四箱与八箱：一行观测仍是一条带 {#statistics-density}

两端按各自完整长度归一化中点后分箱，并结合序列对、方向及 group_by 分组。只将 bins 从 4 改为 8，原始 HSP 行数不变，局部计数的分组粒度改变。bundle_n 是组内观测条数，不是测序覆盖深度；after_stat 读取计算结果。

{{code:comparedensityfour}}
{{figure:comparedensityfour}}
{{code:comparedensityeight}}
{{figure:comparedensityeight}}


## 聚束真正合并显示成员 {#statistics-bundle}

stat_ribbon_bundle 在相同箱和分组内生成汇总带，min_bundle 至少为 2；较小的组仍保留为独立带，不能把全部输入行数直接当作聚束数。显示外包范围不等于重新计算的核苷酸比对，不制造新 identity。bundle_ggchord_ribbons 返回可保留的表和 report，适合先检查、再画图。

{{code:comparebundle}}
{{figure:comparebundle}}
{{code:comparebundlesummary}}
{{output:comparebundlesummary}}


## 统计变量、权重与筛选顺序 {#statistics-parameters}

| 参数/变量 | 默认或含义 |
| --- | --- |
| bins | 80；正整数，每条序列的归一化中点箱数 |
| weight | length；length/count/pident，表示汇总权重来源 |
| group_by | NULL；声明独立分组的列，避免混合不同条件 |
| min_bundle | 2；仅聚束，至少 2 条成员才进行聚束 |
| bundle_n | 分组内原始成员数 |
| bundle_weight | 组内权重总和；count 每行权重为 1，length 用原始长度，pident 用 length × pident / 100 |
| density | 同组权重除以该有序序列对的最大组权重，范围 0–1；不是绝对覆盖深度 |

筛选按阈值选择已有记录；去重选择代表；merge 合并相容邻域；bundle 按归一化箱汇总。先后顺序会改变输入集合，必须保存每步 report 与成员关系。不能把合并后的范围再当原始比对做无说明的生物学推断。



## 用自己的比对与练习 {#statistics-adapt}

把 bins 改为 6，观察相同位置的颜色是否变化，再检查原始行是否都还保留。需要按实验或样本分开时，在源表保留真实条件列并传 group_by。连接统计用于弦图或共线性图；单序列圆图不支持连接层。继续[比对处理](ribbons.html)、[真实案例](real-data.html)和[标度家族](scale-families.html)。
