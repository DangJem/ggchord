## 密集而完整的注释输入 {#complex}

本页独立使用 pETDuet-1，说明完整注释、引物和酶切文字在正常导出尺寸下如何组织；不依赖 pUC19 页的中间对象。

{{code:courseblockg}}

{{code:referencebuilder}}

```r
reference <- make_reference_map(plasmid_example_pETDuet_1)
p <- reference$plot
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogduet}}

## 本图的完整范围 {#source-profile}

20 个特征、5 个独立参考引物和所选集合的 58 条位点记录全部入图，对应 54 个不同切点位置。源文件保存空酶集合，因此明确采用商业非冗余酶中的唯一六碱基及以上识别位点。选择酶集合与任意抽样显示是不同操作。

## 标签与几何的职责 {#layout}

堆叠 Position 安排特征轨道；标签层使用可用内侧空间，必要时外排。真实来源位置保持不动。同一位置的若干酶记录可共享一个标注，不能用文字条数反推来源记录数。

## 增加空间而不删数据 {#size}

同一个图对象可以导出到更宽画布，只改变物理空间，不改变特征、所选位点、链和完整长度。选择尺寸时结合标签诊断查看结果。

```r
ggsave("pETDuet-complete.png", p, width = 16.5, height = 8, units = "in", dpi = 144)
```

## 相关功能 {#next}

[特征堆叠](circular-features.html) · [标签布局](circular-labels.html) · [13 个完整参考图](plasmid-catalog.html)
