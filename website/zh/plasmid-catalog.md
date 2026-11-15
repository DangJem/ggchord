数据来源、展示子集与重建规则见[旗舰数据来源与构建](flagship-data.md)。本页直接加载随包输入，集中设置绘图参数。

这些图保留各质粒的完整来源特征、独立参考引物和所声明酶集合的全部位点。没有等距抽取三个切点。每个图使用自己的完整序列、同一套标注图层和 10 × 8 in 尺寸；点击查看高清图。

<a id="prepare"></a>

## 准备与完整显示规则 {#prepare}

{{code:courseblockg}}

{{code:referencebuilder}}

`reference$features`、`reference$primers`、`reference$sites` 是实际入图的完整表；`reference$all_site_count` 记录原始候选数，完整搜索结果可按构建脚本重建，`reference$policy` 说明使用的酶集合。特征几何不删减，位点集合不做固定数量抽样。

<a id="counts"></a>

## 输入与识别结果 {#counts}
{{plasmidcatalog}}

<a id="catalogpsb"></a>

## pSB1C3：小载体与基本骨架 {#catalogpsb}
{{code:courseblockh}}

{{figure:catalogpsb}}

<a id="catalogpuc"></a>

## pUC19：常见注释与复制起点 {#catalogpuc}
{{code:courseblocki}}

{{figure:catalogpuc}}

<a id="catalogpbr"></a>

## pBR322：多个抗性区与方向 {#catalogpbr}
{{code:courseblockj}}

{{figure:catalogpbr}}

<a id="catalogblue"></a>

## pBluescript II SK(+)：局部调控与短特征 {#catalogblue}
{{code:courseblockk}}

{{figure:catalogblue}}

<a id="catalogpet"></a>

## pET-28a(+)：表达区与标签 {#catalogpet}
{{code:courseblockl}}

{{figure:catalogpet}}

<a id="catalogduet"></a>

## pETDuet-1：双表达单元，进入密集案例 {#catalogduet}
{{code:courseblockm}}

{{figure:catalogduet}}

<a id="catalogdonr"></a>

## pDONR221：重组相关结构 {#catalogdonr}
{{code:courseblockn}}

{{figure:catalogdonr}}

<a id="catalogpcdna"></a>

## pcDNA3.1(+)：表达载体的注释层级 {#catalogpcdna}
{{code:courseblocko}}

{{figure:catalogpcdna}}

<a id="catalogtre"></a>

## pTRE-Tight-BI：双向调控布局 {#catalogtre}
{{code:courseblockp}}

{{figure:catalogtre}}

<a id="catalogtripz"></a>

## pTRIPZ：更长序列的显示密度 {#catalogtripz}
{{code:courseblockq}}

{{figure:catalogtripz}}

<a id="catalogcambia"></a>

## pCAMBIA1300：较大植物载体的结构 {#catalogcambia}
{{code:courseblockr}}

{{figure:catalogcambia}}

<a id="catalogearley"></a>

## pEarleyGate 201：另一植物载体的类别对照 {#catalogearley}
{{code:courseblocks}}

{{figure:catalogearley}}

<a id="catalogcas"></a>

## pSpCas9(BB)-2A-GFP (PX458)：大型表达载体与文字适配 {#catalogcas}
{{code:courseblockt}}

{{figure:catalogcas}}

<a id="reference-selection"></a>

## 酶集合与来源边界 {#reference-selection}

保存的 reference profile 优先，父集合为 `commercial_nonredundant`。pETDuet-1 的源文件保存 None、原参考 PNG 却有位点，因此明确使用 `unique_6plus`，不声称恢复了原 PNG 的导出设置。该集合中的每一行均保留；“完整”指这套明确的集合，不等于把所有酶库中的全部候选同时排到图上。

独立引物取 `set = "reference"`。来源没有独立参考引物时计数为 0，不补造记录。具名子段保存在 `segments` 中，不是额外的顶层特征。相同位置的酶可共用标注，但原始记录仍在表中。

<a id="next"></a>

## 按用途查看 {#next}

[pUC19 完整图层说明](plasmids.md) · [pETDuet-1 密集布局](circular-dense-case.md) · [标签与尺寸](circular-labels.md)

<span id="gallery-limits"></span>
