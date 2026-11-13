## 研究问题与输入身份 {#question}

pUC19 案例先回答如何把一条完整真实 DNA、已识别注释、引物命中与酶切位点分别成层；密集 pETDuet-1 案例再回答在固定导出尺寸下怎样保留可读层级。完整长度来自各自的原始 FASTA，不从最末注释猜测，也不靠改名字把 pUC19c 当作 pUC19。

特征目录的来源与方法、引物目录与精确搜索、限制酶数据库及甲基化背景保留各自元数据。参考目录的整理注释不是未知样本的实验验证。以完整 DNA 精确匹配参考时才能使用对应参考集合，其他输入需检查搜索结果。

## 使用有来源的序列 {#input}
随包 `plasmid_example_pUC19` 原样来自仓库 `examples/plasmid/pUC19.fna`，不会为显示
修改碱基或序列标识。本流程使用包中有来源记录的识别目录，不需要网络请求，也不
为图形伪造实验测量值。

{{code:plasmidinput,frontplasmidsummary}}

{{output:frontplasmidsummary}}

先查看 `plasmid_features`、`plasmid_primers` 与 `display_sites`，区分完整搜索与显示选择。

## 组合独立图层 {#draw}
下面按骨架、注释、标签、引物和酶切的顺序添加图层。注释与标签共享同一个 Position；每个图层只读取自己的数据。

{{code:plasmid}}
{{figure:plasmid}}
{{results:plasmid}}

先看“绘制酶切位点”一行。本例只把 3 个选定切点交给图层。完整搜索结果保存在 `plasmid_sites`，可以回查未画出的位点。

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
本例精确搜索 universal 目录，参考引物注释与精确寡核苷酸搜索有不同来源。
近似蛋白注释识别遵循[识别专题](recognition.html)的限制，不验证生物学功能或
引物实验。文字可以移位、在显式字号下限内缩小或隐藏。应保留几何并在最终尺寸
检查导出标签原因，而不是要求每条注释必须显示文字。

## 应用到自己的数据 {#adapt}
将输入替换成有明确名称的 DNA 序列，先检查搜索结果再添加几何层。已有整理注释时，
直接提供注释表，无需再次运行识别。[堆叠与形状](annotations.html)说明几何，
[标签](labels.html)说明适配，[输出](output.html)说明设备与尺寸。

`geom_primer_label_repel(side = "inside")` 将引物引出标签放在骨架内侧。内侧酶切标签的几何与导出的注释布局表都记录内侧位置。

完整内侧场景仍可能出现长引物文字与其他 feature 区域交叠；默认
`max_overlaps = Inf` 不认证跨 feature 图层的零交叠。空间不足时优先使用默认外侧，
或通过 `show_location = FALSE` 缩短显示文字，并在最终尺寸检查结果。

## 完成后检查 {#check-result}
从中心的名称与长度开始，沿骨架读位置，再读类别色、箭头方向和外侧标签。核对 display_sites 的酶名；删减显示不会改变 plasmid_sites。

## 保存本页结果 {#save-result}
在工作目录中保存刚刚创建的图；预览时使用相同宽高。

{{code:courseblocku}}

## 保存前的证据检查 {#evidence-check}

按图中已声明的宽高调用 view_ggchord，再用相同尺寸 ggsave 保存。核查 ID、源链、比对端点方向、筛选报告、图例单位及标签隐藏原因；PNG 能渲染不代表字体在所有 PDF/SVG 阅读器中一致。保留完整源表、显示子集、参数、软件版本和输入哈希，使图可以追溯到记录。数据来源说明可下载自[资产与来源记录](../assets/ASSET-SOURCES.md)。

## 来源权利与复现 {#source-rights}

包代码采用 MIT 许可；原始数据库记录及参考目录保留各自使用条件，包许可不重新授权第三方数据。导出时应保留记录身份和来源。[资产与来源记录](../assets/ASSET-SOURCES.md)、输入归档哈希及本页准备脚本共同描述本例，计算识别结果仍需实验核实。
