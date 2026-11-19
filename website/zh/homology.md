<a id="membership"></a>

## 读取已有成员 {#membership}

随包文件保存公开植物数据的已有分组。每一行给出组和成员。

```r
groups <- system.file("extdata", "plant-feature-groups.txt", package = "ggchord")
members <- read_feature_relations(groups, format = "named_groups")
head(members$membership[c("group_id", "feature_id")])
```

<a id="links"></a>

## 确认两端的来源 ID {#links}

成对关系表保留拟南芥与大豆的来源 ID。提供带坐标的特征表后，才能为这两端找到区间。

```r
pairs <- system.file("extdata", "plant-feature-pairs.homology", package = "ggchord")
positions <- system.file("extdata", "plant-pair-features.gff", package = "ggchord")
anchors <- read.table(positions, col.names = c("accver", "feature_id", "start", "end"))
joined <- read_feature_relations(pairs, feature_data = anchors)
head(joined$ribbon_data[c("qfeature_id", "sfeature_id", "qstart", "sstart")])
```

<a id="focus"></a>

## 选择局部窗口 {#focus}

为关系画图时，先独立取得完整序列长度，再根据两端的锚点选择窗口。见[局部视图](views.md)。

<a id="projection"></a>

## 蛋白与基因组坐标 {#projection}

已有成员关系标识相关特征。需要精确的氨基酸到核酸坐标时，见[有 CDS 的投影示例](projection.md#cds)。

<a id="check-result"></a>

## 检查关系 {#check-result}

检查源 ID 匹配以及各端的坐标单位。输入第三列的原始分数保留其来源含义。

<a id="save-result"></a>

## 导出已匹配的关系 {#save-result}

返回对象的 ribbon_data 是普通表，可用于另存带两端区间的数据。
