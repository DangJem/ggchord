# 字符与语义细节原型结论

这是v10 §4.9–4.10、§10.13的本版研究交付；不导出新字符Geom、AlignmentView或magnify入口。
代码在semantic-prototypes.R，真实来源及坐标验证在check-semantic-detail.R。

## 两层能力与数据边界

A层：实际字符片段需要accver、完整字符串、明确DNA/RNA/AA字母表及原始位置。
本版通过read_fasta_lengths(keep_sequence=TRUE)提供真实输入，focus同步裁剪；
隔离source_letters生成source_accver/source_position、local/display_position及原始/显示字符。
DNA反向互补保留降序source_position；AA反向仅反转，不作“蛋白互补”。
字母表显式声明，保留大小写；gap不属于无缺口源片段，不把它计成真实bp位置。

B层：对齐片段需要alignment_column与两侧nullable source_position，gap为无源位置，
不能插值制造端点之间的错配。未来保留operation、source row、strand、逐列字符及原始
CIGAR/对齐文本；source→column可能因gap/insertion不是简单一对一。
M、=、X、I、D的消费方向必须随输入格式核对；RNA、A3M小写插入、CDS phase及翻译表
不能通过统一删字符或换大小写掩盖差异。本版不执行比对/翻译推断。

v10列出的ggmsa/ggseqalign等应作为成熟MSA表达与输入生态的比较对象，ggseqlogo为
频率/信息量Logo的可选对接方向，Jalview/ESPript/Geneious为静态阅读与注释层级的设计参考。
这沿用任务书提供的比较目标，不代表本轮重测这些外部工具或复制其实现。
本版选择先补真实字符、坐标和物理空间这一共同缺口，不重建完整MSA/Logo生态或编辑器。

## 实际原型

- 从原始两个核酸FASTA各取101–160的60bp，从四蛋白文件前两条各取21–60的40aa，
  共200个source-position记录。原文件哈希在结果JSON中，CSV与独立静态图保留字符与坐标。
- 图中AA与DNA采用独立原始位置轴，声明这些是独立片段，列位置相同不代表同源对齐。
- 实测等宽字体3mm时，60位点在180mm可用宽度下允许字母，在12mm下要求tiles。
  预算不缩小字号；原型只输出退化决策，不声称自动接入所有ggchord图层。
- 合成字母表用例只验证IUPAC含糊码与大小写的可逆方向；不当作真实生物观测。
- 两次focus将101–160再取局部11–30，得到原始111–130，字符串及源位置一致。

## M0：局部重建而非裁剪图片

真实pETDuet-1全长5420bp，选择101–1900、1800bp。随包字符串与原始FASTA完全一致。
局部保留10条非primer Feature、2条真实引物结合区间、3条真实识别区间；source重叠边界逐条验证。
跨原点的5404→2/3记录与该窗口无交集，在进入线性focus前显式排除，不使用min/max把它扩为长特征。
主图高亮真实弧段；细节图独立Collinear，显示反向骨架以让源首末与弧端的两条leader不交叉。
显示翻转不改变原始链/字符串，轴显式标注101与1900等源位置。

Primer与restriction采用通用Feature图层表示真实结合/识别区间，不冒充圆图专用引物或切口Geom
已经支持Collinear。3个接近的酶识别记录各自保留图形，共享一个组合名称标签。
Feature标签按明确规则选CDS（排除文字ATG）与rep_origin，另显示两引物名和酶名称；
未声称为全注释自动排版或无隐藏的质粒质量认证。原始15条局部注释在annotations.csv可追溯。
没有信号Stat，因此不存在暗中local/global重算；以后加入count/bin/coverage必须独立声明统计口径。

Composer重用Fusion独立panel、实际mm上下文与连续source frame，在当前设备重建两条装饰leader；
不生成伪HSP，导出的biological connections始终为0。最终12×8英寸PNG已人工查看；
10×7只验证物理重建及来源，不作为第二张视觉验收图。未进行PNG/PDF/SVG发布设备矩阵。

## 原型结论与下一门槛

真实字符A层和M0的来源/局部重排在有限条件下可行。当前M0仅一个选区、单一实例、
Circular中心+Collinear右侧、原生identity位置标度；物理投影使用本版的共享导出上下文。
它不是自动inset放置、任意scale/Coord透传或完整Adaptive。独立gTree由ggsave保存，
源图仍为可编辑ggplot对象；研究gTree不应成为稳定公开API。

下一步须先准备两侧真实字符与gapped对齐证据，再验证alignment-column双映射、
负链gap、CDS翻译/phase和codon边界；之后评审共享字符样式与语义放大入口。
M1–M3、多窗口/跨原点/嵌套、任意scale及自由Fusion仍有独立门槛，未纳入本版承诺。
PDF/SVG和完整多尺寸/字体质量矩阵待正式发布授权，不用本地原型替代。
