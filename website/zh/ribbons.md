## 从原始连接表开始 {#input}
先用[真实记录工作流](real-data.html#selection)创建 `seq_real`、`link_real` 和
`real_filtered`，保留原表及筛选报告。

{{code:ribbons}}
{{results:ribbons}}

## 理解每步转换 {#transformations}
筛选保留满足已声明阈值的观测；互惠重叠去重选择已有代表观测，不累加相似度；
合并把邻近兼容区间组合为显示外包范围；聚束按邻近 query/subject 位置形成汇总。
逐步检查 `$report` 和源行成员，特别注意冲突的额外元数据。
合并或聚束范围不是重新计算的核苷酸比对，估计长度保留标记，不凭空生成相似度。

## 密度是显式统计 {#density}
{{figure:density}}

`stat_ribbon_density()` 映射计算后的密度，`stat_ribbon_bundle()` 在 build 时执行
明确聚合，普通 ribbon geom 仍逐行绘制。计数、已报告长度和相似度加权有各自单位，
更不透明的连接不一定表示更强的独立同源证据。没有方向证据的点锚，不能因端点
相等就获得聚合方向。

## 选择区间 ribbon 或点连接 {#points}
先运行[共享准备代码](index.html#quick-start)：

{{code:pointlinks}}
{{figure:pointlinks}}

`support` 是已声明的合成演示值。line 使用 `qpos/spos`，而非区间端点，箭头由
query 指向 subject。line 默认避障为 `none`，ribbon 为 `uniform`。
全局弦图旋转属于坐标系，直线模式选对与相邻行默认见[三种模式](modes.html#collinear)。

点连接的颜色条明确指定 `available_aes = "link_colour"`。普通 ggplot2 颜色条和
包装器默认值不接受全部扩展角色；自定义 guide 时应明确选择对应美学角色。
