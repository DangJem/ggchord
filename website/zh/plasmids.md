数据来源、展示子集与重建规则见[旗舰数据来源与构建](flagship-data.md)。本页直接加载随包输入，集中设置绘图参数。

<a id="input"></a>

## pUC19 的完整输入 {#input}

本页说明一张完整参考图的图层构成。pETDuet-1 的密集标注与尺寸选择见[独立页面](circular-dense-case.md)。

{{code:courseblockg}}

{{code:referencebuilder}}

```r
reference <- make_reference_map(flagship_circular_example$references[["pUC19"]])
p <- reference$plot
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogpuc}}

<a id="layers"></a>

## 阅读各图层 {#layers}

双线骨架提供来源位置；全部参考特征保留原区间、链和分段。特征标签与几何共享堆叠 Position；独立引物和限制位点分别使用自己的数据表与图层。中心文字给出序列名称和完整长度。

<a id="sites"></a>

## 检查完整数据表 {#sites}

```r
reference$policy
nrow(reference$features)
nrow(reference$primers)
nrow(reference$sites)
```

该来源 profile 包含 9 个特征、39 条限制位点记录，对应 33 个不同位置；没有独立参考引物行。原始候选数保留在 `reference$all_site_count` 中，完整搜索结果可按构建脚本重建，入图的参考集合没有再抽取三个位点。

<a id="adapt"></a>

## 使用自己的注释 {#adapt}

已有特征、引物和位点表可以直接交给独立图层，无需先做序列识别。完整长度和来源单位仍须独立提供。

<a id="save-result"></a>

## 按相同尺寸保存 {#save-result}

```r
ggsave("pUC19-complete.png", p, width = 10, height = 8, units = "in", dpi = 144)
```


<span id="check-result"></span><span id="draw"></span><span id="evidence-check"></span><span id="limits"></span><span id="question"></span><span id="source-rights"></span>
