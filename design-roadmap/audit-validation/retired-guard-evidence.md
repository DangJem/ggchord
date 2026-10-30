# 迁移守卫的有限直接验证

基线：`08f753d5b1061fdbaa9ee4eeb91b41f5c8c6782f`，R 4.6.0 / ggplot2 4.0.3。
`audit-reviewed-details/retired-entries.csv` 登记130条旧参数、旧取值或映射事实，
不等于130个公开函数；明确区分“已移除”和offset“仍接受但弃用”。

单一最小脚本检查了三个入口触发，未build图形：

- `ggchord(data.frame(accver="A", length=1000), rotation=0)` 明确提示已移除
  rotation，并给出 `coord_chord(rotation = ...)`。
- `geom_gene_label_repel(gene_label_layout="aligned")` 明确提示已移除取值，
  并给出 `gene_label_layout="auto"`。
- `geom_gene_label_repel(gene_label_segment="solid")` 确实停止，但R先报
  “argument 1 matches multiple formal arguments”，尚未进入函数体中的raw-call守卫。
  不能宣称该旧前缀总能给出本包自定义迁移文字；需要改签名或其他匹配策略时在
  v0.16.1比较方案，不为错误提示重新引入已移除功能。

130条记录未逐条执行；历史删除版本以既有NEWS类别表与明确错误文字为证据，
没有声称已逐tag查出每项最早删除时间。未运行testthat、包检查或发布设备验收。
