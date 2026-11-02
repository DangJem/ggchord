<div class="task-intro"><p>用随包 pUC19 DNA 完成一张可阅读的质粒图：注释、引物和酶切分别成层，显示选择与完整结果分开保存。</p></div>

## 本页准备 {#prepare}
在新的 R 会话中先运行下方代码，再按页面顺序继续。配色和排版是可选的显示设置，输入只选择需要的列。

{{prepare:featurepalette,docstyle}}

其中定义了本例使用的配色。

## 使用有来源的序列 {#input}
随包 `plasmid_example_pUC19` 原样来自仓库 `examples/plasmid/pUC19.fna`，不会为显示
修改碱基或序列标识。本流程使用包中有来源记录的识别目录，不需要网络请求，也不
为图形伪造实验测量值。

{{code:plasmidinput}}

先查看 `plasmid_features`、`plasmid_primers` 与 `display_sites`，区分完整搜索与显示选择。

## 组合独立图层 {#draw}
下面按骨架、注释、标签、引物和酶切的顺序添加图层。注释与标签共享同一个 Position；每个图层只读取自己的数据。

{{code:plasmid}}
{{figure:plasmid}}
{{results:plasmid}}

混合的圈外标签按锚点顺序共同排版，包括少量酶切标签与引物共存的情况。
本例使用质粒专用主题和默认 feature 适配；默认先尝试内部文字，再使用相邻或
外侧位置。不要通过不断缩小文字或强制手动偏移掩盖交叉；请在相同导出尺寸下检查。


原始注释保留在 `plasmid_feature_source`。显示表去掉 `primer_bind` 注释，交由独立
引物层表达，避免重复；注释层只接收几何、类型和注释文字。

## 阅读各注释层 {#layers}
双线骨架建立单条圆形序列，通用注释使用 plasmid 预设和显式 stacking Position，
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
近似蛋白注释识别遵循[识别专题](recognition.html)的限制，不验证生物学功能或
引物实验。文字可以移位、在显式字号下限内缩小或隐藏。应保留几何并在最终尺寸
检查导出标签原因，而不是要求每条注释必须显示文字。

## 应用到自己的数据 {#adapt}
将输入替换成有明确名称的 DNA 序列，先检查搜索结果再添加几何层。已有整理注释时，
直接提供注释表，无需再次运行识别。[堆叠与形状](annotations.html)说明几何，
[标签](labels.html)说明适配，[输出](output.html)说明设备与尺寸。

`geom_primer_label_repel(side = "inside")` 将引物引出标签保留在骨架内侧，
包括同时存在外侧酶切及注释标签时。尚未发布的 0.16.6.9000 checkout 已修复此行为；
v0.16.5 可能将请求内侧的标签移到外侧。外侧默认布局和原始坐标不变。

内侧酶切标签继续保留内侧几何，导出 annotation registry 现在也正确记录内侧；
早期版本将这些诊断行误标为外侧资源。

完整内侧场景仍可能出现长引物文字与其他 feature 区域交叠；默认
`max_overlaps = Inf` 不认证跨 feature 图层的零交叠。空间不足时优先使用默认外侧，
或通过 `show_location = FALSE` 缩短显示文字，并在最终尺寸检查结果。

## 完成后检查 {#check-result}
从中心的名称与长度开始，沿骨架读位置，再读类别色、箭头方向和外侧标签。核对 display_sites 的酶名；删减显示不会改变 plasmid_sites。

<div class="next-step">继续：<a href="labels.html">放置和调整标签 →</a></div>

## 保存本页结果 {#save-result}
在工作目录中保存刚刚创建的图；预览时使用相同宽高。

```r
ggsave("pUC19-map.png", p_plasmid, width = 10, height = 9, units = "in", dpi = 144)
```

## 完整示例：pETDuet-1 密集质粒图 {#complex}

先用上面的 pUC19 掌握图层，再用 pETDuet-1 展示复杂注释。以下来自随包真实序列、常见特征与参考引物表；酶切输入保留 commercial_all 中唯一的六碱基及以上识别位点，不把筛选后的图误称为全部酶的所有切点。

{{code:showcaseplasmiddata}}

{{code:showcaseplasmid}}

{{figure:showcaseplasmid}}

20 个特征、5 个参考引物、58 个合并切点构成输入。默认 feature 配色保留紫红编码区域、浅绿抗性基因、黄色复制起点和蓝色调控注释；引物为紫色，与本站其他图的相同语义保持一致。密集标签对尺寸敏感，下载原图或打开大图阅读。

显式视野将密集图固定在文档给出的导出尺寸。`label_fit = "shrink"` 允许短特征文字缩至 2.5 mm；仍放不下的文字保留为外置标签。预览和保存应使用相同的 16.5 × 8 in 尺寸。

```r
ggsave("pETDuet-1-map.png", p_show_plasmid, width = 16.5, height = 8, units = "in", dpi = 144)
```
