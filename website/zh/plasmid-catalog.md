这 13 张图是可选参考图库，不是需要逐页完成的课程。先完成[基础圆图](circular-start.html#input)、[注释与标签](circular-layout.html#input)，再按自己的载体和拥挤程度检索。

这 13 张图使用各自的仓库 FASTA 参考序列。载体用途用于选例，图中注释来自相应的参考序列。完整来源和长度由随包数据及生成规则保留。

## 一次准备全部案例 {#prepare}
共享脚本使用 dplyr 选择序列列。先运行 `library(dplyr)`。

{{code:courseblockg}}

{{recipe:catalog}}

将脚本放入工作目录并运行，然后选一张图跟做。保留全部特征。名称先在圈内放置，放不下的可移到外圈。外圈单独标注最接近三个等距目标位置的实际酶切位点，显示酶名与位置。

## 输入与识别结果 {#counts}
{{plasmidcatalog}}

## pSB1C3：小载体与基本骨架 {#catalogpsb}
{{code:courseblockh}}

{{figure:catalogpsb}}

## pUC19：常见注释与复制起点 {#catalogpuc}
{{code:courseblocki}}

{{figure:catalogpuc}}

## pBR322：多个抗性区与方向 {#catalogpbr}
{{code:courseblockj}}

{{figure:catalogpbr}}

## pBluescript II SK(+)：局部调控与短特征 {#catalogblue}
{{code:courseblockk}}

{{figure:catalogblue}}

## pET-28a(+)：表达区与标签 {#catalogpet}
{{code:courseblockl}}

{{figure:catalogpet}}

## pETDuet-1：双表达单元，进入密集案例 {#catalogduet}
{{code:courseblockm}}

{{figure:catalogduet}}

## pDONR221：重组相关结构 {#catalogdonr}
{{code:courseblockn}}

{{figure:catalogdonr}}

## pcDNA3.1(+)：表达载体的注释层级 {#catalogpcdna}
{{code:courseblocko}}

{{figure:catalogpcdna}}

## pTRE-Tight-BI：双向调控布局 {#catalogtre}
{{code:courseblockp}}

{{figure:catalogtre}}

## pTRIPZ：更长序列的显示密度 {#catalogtripz}
{{code:courseblockq}}

{{figure:catalogtripz}}

## pCAMBIA1300：较大植物载体的结构 {#catalogcambia}
{{code:courseblockr}}

{{figure:catalogcambia}}

## pEarleyGate 201：另一植物载体的类别对照 {#catalogearley}
{{code:courseblocks}}

{{figure:catalogearley}}

## pSpCas9(BB)-2A-GFP (PX458)：大型表达载体与文字适配 {#catalogcas}
{{code:courseblockt}}

{{figure:catalogcas}}

## 下一步 {#next}
比较不同密度后，进入[综合实战](circular-dense-case.html#complex)，添加独立引物与酶切图层。

## 比较时注意 {#gallery-limits}

每图使用自己的完整 FASTA 和来源注释，不能把一个参考的特征或切点搬到不同 DNA。切点按三个均匀目标选择最近实际位置，是显示选择；不表示其他切点不存在。名称可能内置、外移或隐藏，全部特征几何仍保留。密集图需按标注尺寸查看，缩小截图后读不清不能据此修改源区间。


## 与原参考图比较时的选择 {#reference-selection}

本页固定取3条位点用于简化展示，不能用位点数判断识别完整性。要保留原记录
的酶集合选择，先对完整结果使用
`filter_restriction_sites(sites, set = "reference", parent_set = "commercial_nonredundant")`。
它复用保存profile，不能保证恢复旧PNG的全部导出设置或同名酶库；例如pETDuet-1
原文件保存None而旧PNG显示位点。独立引物需另加`find_primer_bindings(sequence,
set = "reference")`及primer图层。具名子段在`segments`中，不计为另一个顶层特征；
虚线分段也不自动生成额外文字。密集全量标签需要按输出尺寸重新检查。
