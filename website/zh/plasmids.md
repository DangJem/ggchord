## 使用有来源的序列 {#input}
随包 `plasmid_example_pUC19` 原样来自仓库 `examples/plasmid/pUC19.fna`，不会为显示
修改碱基或序列标识。本流程使用包中有来源记录的识别目录，不需要网络请求，也不
为图形伪造实验测量值。

{{code:plasmid}}
{{figure:plasmid}}
{{results:plasmid}}

## 阅读各注释层 {#layers}
双线骨架建立单条圆形序列，通用 feature 使用 plasmid 预设和显式 stacking Position，
标签共享该 Position。引物弧线表达结合范围，酶切短线表达切割/显示位置，中心标签
给出序列名称和长度。裁剪或筛选显示不改变完整序列的酶位点计数。

## 有意选择显示位点 {#sites}
`find_restriction_sites()` 返回全部匹配。本图明确选择 EcoRI、BamHI 与 HindIII，
使局部酶切簇保持可读；`plasmid_sites` 保留完整搜索结果。只有完整输入序列精确匹配参考时，
`reference` 筛选才使用已保存的参考显示集合。其他任务可选 `unique`、`six_plus`
或明确酶名。绘图前检查 `display_sites` 及甲基化/来源列。参考集合是显示选择，
不表示未知样本的生物学状态已经确定。

## 理解目录和标签的限制 {#limits}
`primer_names` 展示目录条目。本例精确搜索 universal 目录，参考引物注释与精确寡核苷酸搜索有不同来源。
近似蛋白 feature 识别遵循[识别专题](recognition.html)的限制，不验证生物学功能或
引物实验。文字可以移位、在显式字号下限内缩小或隐藏。应保留几何并在最终尺寸
检查导出标签原因，而不是要求每条注释必须显示文字。

## 应用到自己的数据 {#adapt}
将输入替换成有明确名称的 DNA 序列，先检查搜索结果再添加 geom。已有整理注释时，
直接提供 feature 表，无需再次运行识别。[堆叠与形状](annotations.html)说明几何，
[标签](labels.html)说明适配，[输出](output.html)说明设备与尺寸。
