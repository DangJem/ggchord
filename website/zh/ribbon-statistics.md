密度着色与聚束不是同一个操作。先保留原始 HSP，再说明箱数、权重、分组与显示汇总的意义。

## 真实输入与最小字段 {#input}

本页独立使用随包的完整序列长度和原始 HSP，不运行新比对。必需 ID 与两端坐标之外，weight=length 需要 length，weight=pident 同时需要 length 与 pident，不能从端点猜这些统计量。

{{code:comparedensitydata}}


## 四箱与八箱：一行观测仍是一条带 {#density}

两端按各自完整长度归一化中点后分箱，并结合序列对、方向及 group_by 分组。只将 bins 从 4 改为 8，原始 HSP 行数不变，局部计数的分组粒度改变。bundle_n 是组内观测条数，不是测序覆盖深度；after_stat 读取计算结果。

{{code:comparedensityfour}}
{{figure:comparedensityfour}}
{{code:comparedensityeight}}
{{figure:comparedensityeight}}


## 聚束真正合并显示成员 {#bundle}

stat_ribbon_bundle 在相同箱和分组内生成汇总带，min_bundle 至少为 2；较小的组仍保留为独立带，不能把全部输入行数直接当作聚束数。显示外包范围不等于重新计算的核苷酸比对，不制造新 identity。bundle_ggchord_ribbons 返回可保留的表和 report，适合先检查、再画图。

{{code:comparebundle}}
{{figure:comparebundle}}
{{code:comparebundlesummary}}
{{output:comparebundlesummary}}


## 统计变量、权重与筛选顺序 {#parameters}

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



## 用自己的比对与练习 {#adapt}

把 bins 改为 6，观察相同位置的颜色是否变化，再检查原始行是否都还保留。需要按实验或样本分开时，在源表保留真实条件列并传 group_by。连接统计用于弦图或共线性图；单序列圆图不支持连接层。继续[比对处理](ribbons.html)、[真实案例](real-data.html)和[标度家族](scale-families.html)。
