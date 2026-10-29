# API 变化与迁移记录

**状态：v0.16.0 开发中，首批修复已实现并通过局部直接验证；尚未正式发布。**

本文件只列实际变化与有明确标记的待决策事项。没有新增、删除或重命名公开函数，
也没有改变参数默认值；完整发布验收仍须用户明确授权后执行。

## 已实现变化

| 问题 | 旧行为 | 当前行为 | 用户影响与替代写法 | 删除版本 |
| --- | --- | --- | --- | --- |
| F01 | clean 的非 NULL 错误类型可选表被静默丢为 NULL | 明确报出 ribbon_data/gene_data 类型错误 | 传 data.frame/tibble；有意不提供时用 NULL；支持的区间对象先显式适配 | 不适用，未删除参数 |
| F02 | 反向 ribbon 排序报告的 original_value/new_value 相同 | 分别记录交换前后值 | 调用不变；读取 report 时得到真实变更；数据排序和 direction 语义不变 | 不适用 |
| F11 关联修复 | plasmid wrapper 的 shape 单参被误判为双形状参数 | 单参别名正常转发 | 既有 shape 调用恢复；显式两参仍报原冲突 | 不适用 |

### F01：明确输入，而不是丢弃数据

以下示例的 seqs/ribbons 是完整教学输入；不代表实际科研数据。

```r
library(ggchord)
seqs <- data.frame(accver = c("q", "s"), length = c(20, 20))
ribbons <- data.frame(qaccver = "q", saccver = "s",
  qstart = 2, qend = 8, sstart = 15, send = 3)

# 旧错误调用曾静默得到 ribbon_data = NULL；当前明确报错。
# clean_ggchord_data(seqs, ribbon_data = "not a table")

# 显式传入数据；可选表有意缺席时才使用 NULL。
cleaned <- clean_ggchord_data(seqs, ribbon_data = ribbons, gene_data = NULL)
```

`GRanges`、`GInteractions` 等不是 clean 的新增输入类型。需要时先使用对应的
`as_ggchord_feature()` / `as_ggchord_link()` 转成表，再调用 clean；适配器自身的
坐标、链和可选依赖要求仍适用。不要为了绕过报错把本应保留的数据改成 NULL。

### F02：调用不变，修正报告

```r
sorted <- cleaned$report[cleaned$report$action == "sort", ]
stopifnot(nrow(sorted) == 1L,
  sorted$original_value == "qstart=2, qend=8, sstart=15, send=3",
  sorted$new_value == "qstart=2, qend=8, sstart=3, send=15",
  cleaned$ribbon_data$direction == "reverse",
  ribbons$sstart == 15, ribbons$send == 3)
```

修复限定于 ribbon 排序报告；gene 分支原本已能区分交换前后值。默认排序政策与
prepare 的方向保留差异没有在本轮调整，详见 F03 待决策事项。

### F11：恢复已有 shape 简写

```r
features <- data.frame(accver = "q", start = 2, end = 8,
  strand = "+", type = "CDS")
p <- ggchord(seqs[1, ]) + geom_seq() +
  scale_feature_fill_manual(values = c(CDS = "black")) +
  geom_feature_plasmid(data = features, shape = "block") + coord_circular()
built <- ggplot2::ggplot_build(p)
stopifnot(all(built$data[[2]]$feature_fill == "black"))
```

同样可以继续使用 `feature_shape = "block"`；本轮没有决定哪个名称长期保留。
`geom_feature_plasmid(feature_shape = "arrow", shape = "block")` 仍明确报错，
没有采用静默优先规则。默认几何、布局与用户 scale 优先行为保留。

## 尚未实施的迁移

F03 默认差异、F07 getter/export、F09 gene/feature、F10 既有偏移弃用、F11 唯一
形状入口、标签及连接命名、颜色别名和重导出均尚未新增删除决定。后续若改变，
必须补入旧/新调用、替代方案、影响、错误提示、用户确认状态及确切删除版本；
不能将待确认的删除日期写成承诺。详细方案见 [决策台账](API-AUDIT-DECISIONS.md)。
