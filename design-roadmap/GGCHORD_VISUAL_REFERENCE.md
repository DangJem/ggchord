我挑选了13个最具有代表性的质粒序列，放到了 `examples/plasmid`目录下，每个质粒分别有其在SnapGene官方提供的图谱截图、官方提供的SnapGene软件的 `.dna` 工程文件以及fasta数据。

这里有几个很重要的结论先说在前面：

* 这 13 个文件合计包含 **201 个 feature、218 个 feature segment**，覆盖 **14 种 SnapGene feature type**；另外 **pSB1C3 和 pETDuet-1 还有真正的 Primer 注释**，它们并不属于 Feature 表。
* 这批数据非常适合作为 ggchord 的视觉基准，因为它不仅覆盖 feature 类型，还覆盖了**稀疏/拥挤、长/短、单向/双向、跨圆周零点、重复 cassette、外置 feature 标签、极密 restriction labels、primer labels** 等情况。
* 从原始 `.dna` 可以直接得到 feature 的精确颜色，因此下面的 HEX 色值不是从截图“吸色”猜出来的。
* **不要把这些颜色理解成“SnapGene 中某一种 feature type 永远只能有一种颜色”**。尤其 `CDS` 和 `misc_feature` 明显采用语义化颜色。同一种 `CDS` 可以是浅绿、酒红、粉紫、亮绿、红色或橄榄色。
* 对 ggchord 来说，目标应该是**视觉行为等效**，而不是要求 feature/site 数量与 SnapGene 官方数据库逐条一致。你的数据库多几个或少几个都不重要；同样的数据进入绘图层以后，其**颜色、层次、几何形态、标签位置、避让和连接线行为**应尽可能接近 SnapGene。

下面这份可以直接整理给 Codex。

---

# ggchord：基于 13 个参考质粒总结 SnapGene 圆图视觉规则

## 一、总体目标

以以下 13 个 SnapGene 参考质粒作为视觉 benchmark：

`pBR322`、`pUC19`、`pBluescript II SK(+)`、`pSB1C3`、`pET-28a(+)`、`pETDuet-1`、`pcDNA3.1(+)`、`pTRE-Tight-BI`、`pSpCas9(BB)-2A-GFP (PX458)`、`pDONR221`、`pCAMBIA1300`、`pEarleyGate 201`、`pTRIPZ`。

目标不是逐一硬编码这些质粒，而是从这些实例中归纳出一套**通用 SnapGene 风格圆图语法**。

绘图背景统一使用：

**`#FFFCF5`**

注意：原始 `.dna` 中定义为白色的 feature 应继续保持真正的 **`#FFFFFF`**，不要因为背景改成 `#FFFCF5` 而同步改成米白。这样白色 promoter / terminator 等 feature 才会在背景上保留非常轻微的视觉层次。

---

# 二、SnapGene 总体布局特征

## 2.1 圆形序列本体

SnapGene 的圆图不是“一条圆线 + 一圈箭头”。

参考图共同具有以下结构：

* DNA backbone 是非常清晰的**双线圆环**；
* 主轮廓接近深灰黑色，截图中主要深色接近 `#252525`；
* 两条圆线之间保持很窄的间隔；
* feature 基本位于 DNA 圆环**内侧**；
* restriction-site labels 基本位于圆环**外侧**；
* 需要外置的 feature label 也进入圆环外侧；
* 中心只放质粒名称和长度，不塞入 legend；
* 整体非常强调白/空白区域，不追求把圆内部填满。

质粒名称：

* 居中；
* 粗体；
* 黑色；
* 下一行为普通字重的 `xxxx bp`；
* 名称与长度之间行距很小。

---

## 2.2 圆图大小不是固定的

这一点非常重要。

SnapGene 并不是：

> 每张图固定一个半径，然后强行把所有文字塞进去。

从 13 张图能明显看到：

* `pBR322`、`PX458` 等外部标签较少时，圆可以很大；
* `pETDuet-1` restriction labels 极端密集，因此序列圆明显缩小，给外围标签释放空间；
* 当一侧标签特别多时，圆心也可以视觉上偏向另一侧；
* 最终 canvas 尺寸跟内容复杂度有关。

也就是说：

> **map radius、map center、外围留白应该由 annotation density 共同决定。**

不应为了保持固定半径而允许文字碰撞。

SnapGene 自己的图像生成接口也明确说明，实际输出尺寸可能为了确保文字元素完整显示而发生调整。([SnapGene支持][1])

---

## 2.3 坐标刻度是自适应的

13 个质粒已经表现出不同间隔：

* ~2 kb 质粒可以采用约 250 bp；
* ~3–5 kb 常见 500 或 1000 bp；
* ~9 kb 常见 1000 bp；
* ~12–13 kb 可以达到 2000 bp。

所以不要固定为“每 1000 bp 一个刻度”。

原则应该是：

> 根据 sequence length 和圆周可用空间选择整洁的 nice interval，使一圈大约维持有限数量、可读的主刻度。

刻度特点：

* 刻度从圆环向内；
* 数字也位于圆环内侧；
* 数字沿圆周切向排列；
* 0/起点在约 12 点方向用一条明显更长、更粗的径向标记表示；
* 通常不额外写一个普通的 `0` 刻度文字。

---

# 三、这 13 个质粒分别应该让 ggchord 学什么

| 质粒                       | 主要学习目标                                                                                         |
| ------------------------ | ---------------------------------------------------------------------------------------------- |
| **pBR322**               | 最低复杂度 baseline；长 CDS、ori、promoter、misc feature；基础箭头比例、文字贴合、稀疏布局                                |
| **pUC19**                | MCS 周围极密 restriction labels；同位点多酶合并；lacZα 与小 feature 的径向层次                                     |
| **pBluescript II SK(+)** | MCS、primer_bind、T3/T7 promoter、lac operator、lacZα 同区高度嵌套；小 feature layering                    |
| **pSB1C3**               | BioBrick prefix/suffix、terminator；跨 0 点邻域；真正的 Primer 外置标签                                      |
| **pET-28a(+)**           | 6xHis、T7 tag、thrombin、ATG、RBS、operator 等极短 feature；大量 feature label 外置并与 restriction labels 混排 |
| **pETDuet-1**            | 两套重复表达 cassette；跨 0 点 promoter；极端 restriction-site 密度；Primer 与 site 混排                         |
| **pcDNA3.1(+)**          | 中等复杂度哺乳动物表达图；promoter/enhancer/polyA/origin/selection marker 的标准组合                             |
| **pTRE-Tight-BI**        | 双向 promoter；7 个 tet operator；两个 MCS；双向结构与对称布局                                                  |
| **PX458**                | 超长 Cas9 与极短 NLS/T2A/FLAG/ATG/gRNA 同图；极端长度比例；reporter 颜色                                        |
| **pDONR221**             | Gateway attP1/attP2 + ccdB/CmR cassette；protein_bind 的方向性形态                                    |
| **pCAMBIA1300**          | plant/T-DNA；跨 0 点的 MCS 和 lacZα；植物、细菌 backbone 混合                                               |
| **pEarleyGate 201**      | Plant + Gateway；feature 类型/颜色高度丰富；大量短 feature 外置标签                                             |
| **pTRIPZ**               | 最终压力测试；viral + inducible + reporter + shRNA + IRES + WPRE + LTR + operator；外部 site/feature 混排  |

其中可以把：

**pBR322 → pUC19 → pET-28a → pTRE-Tight-BI → PX458 → pEarleyGate 201 → pTRIPZ**

看成由简单到极复杂的核心 visual regression 骨架。

---

# 四、这 13 个 `.dna` 实际包含哪些 Feature Type

解析得到一共 **14 类 feature type**：

| Feature Type    | 数量 | 主要视觉颜色    |
| --------------- | -: | --------- |
| `CDS`           | 65 | 多种语义颜色    |
| `promoter`      | 41 | `#FFFFFF` |
| `misc_feature`  | 25 | 多种        |
| `protein_bind`  | 25 | `#31849B` |
| `rep_origin`    | 23 | `#FFFF00` |
| `primer_bind`   | 10 | `#A020F0` |
| `terminator`    | 10 | `#FFFFFF` |
| `polyA_signal`  |  9 | `#A6ACB3` |
| `RBS`           |  3 | `#A6ACB3` |
| `LTR`           |  2 | `#FFE4C4` |
| `enhancer`      |  2 | `#FFFFFF` |
| `misc_RNA`      |  1 | `#00CCFF` |
| `repeat_region` |  1 | `#FFE4C4` |
| `intron`        |  1 | `noColor` |

这里尤其不能把：

```text
CDS -> 一个颜色
misc_feature -> 一个颜色
```

处理得过于简单。

SnapGene 这批官方质粒明显带有**语义配色**。

---

# 五、从这 13 个文件得到的实际配色体系

## 5.1 CDS

### 抗性/选择标记：浅绿色

**`#CCFFCC`**

包括：

* AmpR
* KanR
* CmR
* TcR
* HygR
* NeoR/KanR
* PuroR
* BleoR
* BlpR

这是这套 reference 中最稳定的语义色之一。

---

### 一般重要蛋白/CDS：酒红色

**`#993366`**

包括：

* lacZα
* lacI
* rop
* Cas9
* ccdB
* rtTA3
* pVS1 RepA
* pVS1 StaA
* SV40 NLS
* nucleoplasmin NLS
* ATG

这是另一个非常强的 SnapGene 视觉主题色。

---

### 小蛋白标签、peptide/tag：粉紫色

**`#CC99B2`**

包括：

* 6xHis
* T7 tag
* S-Tag
* HA
* 3xFLAG
* T2A
* thrombin site

这是 pET-28a、PX458、pEarleyGate 特别重要的一类。

---

### Reporter

EGFP：

**`#05FD14`**

TurboRFP：

**`#FF0000`**

SnapGene 对 reporter 使用非常鲜艳、接近 reporter 名称本身语义的颜色。

---

### 特殊状态 CDS

pTRIPZ 的 `HygR (inactivated)`：

**`#808000`**

因此不要假设所有 resistance feature 都必须浅绿色；annotation 本身可以携带覆盖默认 palette 的颜色。

---

# 六、其他 Feature Type 的颜色

### replication origin

`rep_origin`

**`#FFFF00`**

包括：

* ori
* f1 ori
* SV40 ori
* pVS1 oriV

这是这套数据中最稳定的一种映射。

---

### promoter

**`#FFFFFF`**

41 个 promoter 全部如此。

包括：

* AmpR promoter
* T7 promoter
* T3 promoter
* lac promoter
* lacI promoter
* CMV promoter
* SV40 promoter
* U6 promoter
* UbC promoter
* CaMV 35S promoter
* minimal CMV promoter
* bidirectional TRE promoter
* MAS promoter
* cat promoter 等。

虽然是纯白填充，但有深色描边，所以在 `#FFFCF5` 背景上仍然可以辨认。

---

### terminator

**`#FFFFFF`**

包括：

* T7 terminator
* rrnB T1/T2
* lambda t0
* his operon terminator
* bacterial terminator
* MAS terminator
* OCS terminator
* polIII terminator

---

### enhancer

**`#FFFFFF`**

CMV enhancer。

---

### protein binding / recombination / operator

`protein_bind`

**`#31849B`**

包括：

* lac operator
* tet operator
* attP1
* attP2
* attR1
* attR2
* loxP511

这也是高度稳定的类别颜色。

---

### primer_bind feature

**`#A020F0`**

包括：

* M13 fwd
* M13 rev
* KS primer
* SK primer

注意：

> 这是 **Feature Type = primer_bind**，不等于 SnapGene 独立的 Primer annotation。

后者见后文。

---

### RBS

**`#A6ACB3`**

---

### poly(A) signal

**`#A6ACB3`**

包括：

* SV40 poly(A)
* bGH poly(A)
* CaMV poly(A)

---

### LTR

**`#FFE4C4`**

5′ LTR、3′ LTR。

---

### repeat_region

AAV2 ITR：

**`#FFE4C4`**

---

### misc_RNA

gRNA scaffold：

**`#00CCFF`**

---

### misc_feature

这一类特别不能统一颜色。

MCS / BioBrick prefix/suffix：

**`#99CCFF`**

IRES / TMV Ω / bom：

**`#A6ACB3`**

HIV-1 ψ / RRE / cPPT-CTS / LB-RB T-DNA repeat：

**`#FFE4C4`**

shRNAmir insertion site：

**`#CC99FF`**

WPRE：

**`#FFFFFF`**

因此 `misc_feature` 应当被视为：

> “需要进一步按 semantic identity 决定颜色”的兜底类型。

---

# 七、Feature 的径向层次和重叠规则

这是目前 ggchord 最应该重点学习的部分之一。

从 pUC19、pBluescript、pCAMBIA1300 等图可以非常明显地看到：

> **发生坐标重叠时，长 feature 通常更靠近 DNA backbone，即更靠外；短 feature 依次向圆心方向叠放。**

例如 pBluescript：

* 大的 `lacZα` 靠近圆环；
* MCS、primer、operator、promoter 等短 feature 被排在其内侧。

pCAMBIA1300 也是如此：

* 跨零点的大 `lacZα` 在较外轨；
* MCS、primer/operator 等在里面。

SnapGene 当前官方说明也明确描述其默认行为是：

> overlapping features 按“最长在上、最短在下”的顺序 tile；用户可显式将某些 feature 设成优先显示以打破默认顺序。([SnapGene支持][2])

映射到圆图，应理解成：

**长 feature 优先占靠近序列圆的外层轨道，短 feature 向圆心内缩。**

因此不要简单按：

* 输入顺序；
* feature 类型；
* 正负链；
* 随机 greedy 顺序

分配径向层。

### 应遵守的视觉规则

当两个 feature 不重叠：

> 可以使用同一径向轨道。

当发生重叠：

> 需要增加内层轨道。

默认排序：

> 更长 → 更外；更短 → 更内。

但如果存在“显式优先显示”语义，则它可以突破长度规则。

这 13 个文件没有提供足够证据去定义“长度完全相同时如何 tie-break”，不要为了这些示例虚构特殊规则。

---

# 八、Feature 箭头的总体几何语法

## 8.1 DNA 长度只决定圆周方向长度

Feature 的 bp 长度主要映射为：

> **angular span**

而不是：

> bp 越长，feature 越粗。

Feature 的径向厚度总体保持稳定。

所以 Cas9 即使几千 bp，也只是圆周上占据更长的角度，并不会画成超粗箭头。

---

## 8.2 directional feature

例如 CDS、带方向的 origin、某些 protein_bind。

视觉上是：

**annular arrow / curved ribbon arrow**

特点：

* body 沿圆弧延伸；
* body 的径向厚度基本恒定；
* 3′ 端形成明显三角/chevron arrowhead；
* arrowhead 比 body 本身略宽；
* 有细深灰轮廓；
* fill 使用 feature color；
* 正向沿坐标递增方向；
* 反向沿坐标递减方向。

`.dna` 中：

* `directionality=1` 对应 forward；
* `directionality=2` 对应 reverse；
* `directionality=3` 可表示双向。

pTRE-Tight-BI 的：

**bidirectional TRE promoter**

就是非常重要的双向参考。

---

## 8.3 短 directional feature

不能单纯将大箭头同比缩小。

例如：

* HA
* 6xHis
* ATG
* NLS
* T2A
* attR
* primer_bind

SnapGene 会保证一个**最低可见 glyph 尺寸**。

因此很短的 feature 看起来往往更像：

* 小箭头；
* 小 chevron；
* 小 wedge；

而不是薄得看不见的一个三角像素。

---

## 8.4 non-directional feature

例如：

* MCS
* poly(A)
* terminator
* IRES
* T-DNA repeat
* tet operator
* RBS

一般呈现成：

**弧形矩形 / annular block**

即：

* 无箭头头部；
* 两端大体平直；
* 顺着序列圆曲率；
* 保持相同 radial thickness；
* 仍然有深色细描边。

---

# 九、不同 feature 类型应体现不同 glyph 语义

不应仅仅“所有 feature 都画成箭头”。

参考图体现的倾向是：

### CDS

明显实心 directional arrow。

### promoter

白色 hollow/white arrow。

即使 fill 是 `#FFFFFF`，通过黑/灰轮廓仍表现出方向。

### bidirectional promoter

双向箭头语义。

### origin

有方向时为醒目的黄色 curved arrow。

无方向的 origin 可以呈黄色 block。

### protein_bind

短 cyan/teal block 或 directional wedge。

例如：

* tet operator：短矩形；
* attP/attR：方向性更明显。

### terminator / poly(A) / enhancer

一般采用比较中性的白色/灰色 block/ribbon。

### peptide/tag

粉紫色、小型 arrow/wedge。

### MCS

浅蓝色 block，不应该伪装成 CDS 箭头。

---

# 十、Feature 上的“虚线”是什么

这一点从 `.dna` 数据可以明确确定。

例如 pBR322 的 AmpR：

```text
3293–4084
4085–4153 signal sequence
cleavageArrows = 4084
```

pUC19、pBluescript、pETDuet、pcDNA3.1、pTRE、pTRIPZ 的 AmpR 都具有类似信息。

另外还有：

* pET-28a thrombin site；
* PX458 3xFLAG；
* PX458 T2A。

所以参考图上横穿 feature 的虚线/点状切痕不是：

> “feature 被拆成两个 segment 的通用边界线”

而主要是：

> **cleavage site marker**

即 `.dna` 中的 `cleavageArrows`。

应在 feature 内的具体 cleavage coordinate 上绘制一条：

* 横跨箭头径向宽度；
* 与 feature 轴线近似垂直；
* 点状/虚线或视觉上的切口；

来表示切割位置。

特别注意：

> 多 segment feature 不等于一定画虚线。

例如 lac promoter 可以由 `-35 / core / -10` 多个 segment 组成，但没有对应 cleavage metadata，因此不能给每一个 segment boundary 都加一条虚线。

---

# 十一、多 segment feature

SnapGene 会把属于同一个 feature 的相邻 segments 视为一个整体。

典型例子：

AmpR：

```text
main CDS
+
signal sequence
```

lac promoter：

```text
-35
core
-10
```

EGFP 也存在多个相邻 segment。

总体效果应是：

* 仍然作为一个 feature；
* 连续布局；
* 不人为变成多个互相独立的 feature；
* 可以保留 segment-specific boundary / cleavage 等语义。

---

# 十二、跨越 circular origin 的 feature

这套 benchmark 有很多很好的案例：

### pCAMBIA1300

MCS：

`8959..56`

lacZα：

`8945..219`

### pETDuet-1

T7 promoter：

`5404..2`

### pTRE-Tight-BI

bidirectional TRE promoter：

`2787..318`

这些都必须作为：

> 一个连续跨越 12 点/0 bp seam 的 feature

来处理。

不能出现：

* 末端一截；
* 起点一截；
* 被误判成两个无关 feature；
* 两个独立 label。

这是 ggchord 圆图必须通过的测试。

---

# 十三、Feature label 的文字布局

这是 SnapGene 风格最明显的特征之一。

## 13.1 足够长的 feature：文字沿 feature 圆弧

例如：

* AmpR
* Cas9
* lacI
* ori
* pVS1 RepA
* KanR
* PuroR

标签的视觉方向不是简单地：

> 在 feature 中心放一个旋转后的直线字符串。

而是呈现为：

> **沿 feature 中轴圆弧排布、字形局部方向跟随圆弧切线。**

视觉上文字和箭头共同“弯”在圆周上。

因此 ggchord 应以：

> **arc-following label**

作为目标，而不是简单旋转一次的普通直线文本。

---

## 13.2 标签始终尽量保持可读

对于圆的不同半区：

* 不要允许文字整体倒置；
* 根据所在半圆以及 feature 方向选择合适文字路径方向；
* 必要时翻转 text path；
* 箭头方向仍保持生物学方向，但标签应优先保证人能够正向阅读。

---

## 13.3 字色自动根据 feature fill 对比

参考图中：

浅色 feature：

* AmpR
* ori
* KanR

通常使用：

**黑色文字。**

深色 feature：

* lacZα
* Cas9
* rtTA3
* pVS1 RepA
* ccdB

通常使用：

**白色文字。**

因此文字颜色不应固定为黑色。

---

# 十四、什么时候 feature label 不放在 feature 上

SnapGene 有至少三级策略：

### 第一优先

能放入 feature：

> 放在 feature 内或紧贴 feature，沿弧排列。

### 第二优先

feature 太短，但附近圆内部仍有空间：

> 将文字放在 feature 附近的圆内区域，并保持清楚的空间关联。

### 第三优先

如果内部太拥挤，或者 feature 极短：

> **把 feature label 移到圆外。**

这一点在你的参考图里非常明显。

例如 pET-28a：

* MCS
* thrombin site
* 6xHis
* ATG
* RBS
* lac operator
* T7 promoter

很多都被放到外面。

pTRIPZ：

* HIV-1 ψ
* tet operator × 多个
* minimal CMV promoter
* TurboRFP

也出现大量圆外 feature callout。

pEarleyGate：

* LB T-DNA repeat
* TMV Ω
* ATG
* HA
* attR1
* lac UV5 promoter

同样如此。

SnapGene 官方现在也明确说明：当 circular map 中 feature label 拥挤时，会将部分 feature 名称放到圆外；如果强制禁止圆外标签，则部分标签可能被隐藏。([SnapGene支持][3])

---

# 十五、圆外 feature label 的外观

这是需要 ggchord 特别模仿的一套样式。

它不是普通裸文字。

通常表现为：

**rounded callout**

即：

* 小圆角矩形；
* 黑色文字；
* 填充颜色取 feature 原色的明显浅化版本；
* 边框保持同一 hue；
* 与 feature 通过细灰色 leader line 相连。

例如：

### tet operator

原 feature：

`#31849B`

外置 label 呈明显的浅青蓝色 capsule。

### TurboRFP

原 feature：

`#FF0000`

外置 label 呈浅红/粉红背景，而不是直接用全红底。

### HA / ATG / peptide

原色：

`#CC99B2` / `#993366`

外置 label 使用相应的浅粉色视觉语言。

### promoter

原 feature：

`#FFFFFF`

外置 label 基本就是白色圆角框。

因此：

> 外置 label 应保持 feature identity，但不应该把原始高饱和度 fill 原封不动用在标签背景上。

应该明显 lightening / pastel 化。

---

# 十六、restriction site 标签的基本外观

Restriction sites 与 feature labels 是另一套视觉语言。

特点：

* 位于 DNA 圆外；
* 文字基本水平；
* 黑色；
* 没有彩色 box；
* 通过细灰色 leader line 指回精确 cut coordinate；
* enzyme name 与 coordinate 组合显示。

---

## 16.1 左右两侧的文字顺序不同

这是很容易忽略、但非常 SnapGene 的细节。

右侧：

```text
EcoRI  (952)
```

即：

**enzyme → coordinate**

因为 enzyme name 更靠近圆。

左侧：

```text
(2503)  SspI
```

即：

**coordinate → enzyme**

同样保证 enzyme name 更靠近圆。

所以原则是：

> **enzyme name 始终处于靠近 DNA 圆的一侧；coordinate 位于更外围。**

这比简单固定 `"enzyme (position)"` 更接近 SnapGene。

---

# 十七、相同坐标的多个 restriction enzymes

SnapGene 不会为完全相同的 cut coordinate 拉出五条平行线。

而是合并：

```text
AvaI - BsoBI - KpnI - TspMI - XmaI   (412)
```

或者左侧：

```text
(1300)   BtgI - NcoI - StyI
```

特点：

* 共用一个 coordinate；
* 共用一条 leader；
* enzyme name 用 `" - "` 拼接；
* 整组视为一个 label box 进行 collision avoidance。

这点 pUC19、pBluescript、pET-28a 等特别明显。

---

# 十八、unique cutter 的字体

你的截图中粗体与普通字体不是随机的。

SnapGene 默认：

> **只切一次的 unique cutter 使用粗体。**

官方说明也明确如此。([SnapGene支持][4])

因此：

* unique cutter → **bold**
* repeated cutter → regular

例如 pTRIPZ 中多次出现的 `AanI` 就是普通字重，而很多只有一次的酶是粗体。

---

# 十九、灰色 restriction enzymes 和 `*`

你的截图还能看到：

```text
BsaBI *
BclI *
PflMI *
```

呈灰色。

这里反映的是 methylation 语义。

SnapGene 对被 methylation 阻断的 restriction site 使用灰色显示；历史版本的 release notes 也明确提到“restriction sites blocked by methylation”应显示为 gray。([SnapGene][5])

新版对 methylation-sensitive enzyme 也使用 `*` 提示。([SnapGene][6])

所以 ggchord 如果数据库能提供这类信息，应保持：

* 灰色 enzyme name；
* `*` 标志；
* 但不需要因为数据库差异而强求与 SnapGene 数量完全一致。

---

# 二十、restriction leader line 的几何形态

连接线不是曲线，也不是大量 spline。

总体是：

> **细、浅灰、折线式 polyline。**

参考色接近：

**`#7F7F7F`**

表现形式根据拥挤程度变化。

---

## 20.1 稀疏位置

如 pBR322 的孤立位点：

* 从序列圆精确 cut coordinate 出发；
* 基本径向/斜向外伸；
* 可以接近一条直线；
* 到文字前停止；
* 无箭头头部。

---

## 20.2 密集位置

如：

* pUC19 MCS
* pBluescript MCS
* pET-28a MCS
* pETDuet MCS

采用明显的 fan-out。

视觉结构大致是：

```text
sequence anchor
     │
     │ short radial section
     └──────── oblique leader ───── label
```

即：

1. 从真实 cut coordinate 出发；
2. 先走一小段基本径向线；
3. 再产生 elbow；
4. 向不同高度/方向扇开；
5. 最后接水平 label。

大量近邻位点会形成一个非常整齐的“扇骨”结构。

---

# 二十一、密集 restriction labels 的关键规则

不要仅仅采用局部 repel。

SnapGene 的效果更像是：

> **先依据圆周 anchor 顺序建立整体次序，再分配外围 label 行。**

因此应做到：

* 邻近 cut sites 在文字区域仍保持原来的圆周顺序；
* 尽量不交换上下顺序；
* 因此 leader lines 基本不交叉；
* label 行之间维持近似一致间距；
* 高密度时整体向外扩展；
* 必要时缩小圆，而不是让几十个标签互相覆盖。

这是 pETDuet-1 最重要的 benchmark。

---

# 二十二、Feature 外置 label 与 restriction labels 不是两套独立布局

这是你第 6 点里最关键的问题。

当前：

* feature label → `geom_feature_label_repel`
* restriction site → `geom_restriction_site`

如果这两者各自独立排版，然后最后叠起来，即使两个函数内部各自都“没有重叠”，最终仍然会出现：

> feature labels 与 restriction labels 相互覆盖。

而 SnapGene 明显不是这样。

看 pET-28a 右侧的顺序：

restriction enzyme
→ MCS feature label
→ restriction enzyme
→ thrombin / 6xHis
→ restriction enzyme
→ ATG
→ RBS
→ restriction enzyme
→ lac operator
→ T7 promoter
→ restriction enzyme

它们按照各自在序列上的 anchor **混在同一个外围 label 系统中**。

pTRIPZ 右侧同样非常明显：

restriction sites
→ tet operator callouts
→ minimal CMV promoter
→ restriction sites
→ TurboRFP
→ restriction sites

所以：

> **SnapGene 的外部 feature labels 和 restriction labels 在视觉上属于同一个 global external annotation layout，而不是两个互不知情的图层。**

---

# 二十三、建议 ggchord 采用统一的“外围 annotation layout”

这里不是要求改变用户层面的两个 geom。

用户仍然可以继续调用：

* `geom_feature_label_repel`
* `geom_restriction_site`

但最终视觉排版阶段不能各排各的。

应该把需要放到圆外的 annotation 统一看成一组：

```text
External annotations
├── restriction-site label
├── external feature label
└── primer label（如果显示）
```

每个对象都拥有：

* anchor coordinate；
* anchor angle；
* preferred side；
* bounding box；
* visual style；
* leader-line style；
* priority；
* 是否允许向外移动。

然后**共同计算最终位置**。

---

# 二十四、统一外围布局的推荐流程

### 第一步：Feature label 先判断能否留在圆内

长且空间充足：

> 留在 feature 上。

短但附近有空位：

> 留在圆内、靠近 feature。

明显放不下：

> 转换成 external feature callout candidate。

---

### 第二步：Restriction site 全部生成 external candidates

每一个 site 或同坐标 site group：

* 精确 anchor；
* label text；
* coordinate；
* unique/nonunique；
* methylation state。

---

### 第三步：将所有圆外 annotation 合并

不再区分：

> “这是 feature 图层的位置空间”
> “这是 restriction 图层的位置空间”

它们必须彼此可见。

---

### 第四步：按圆周位置分区

例如左、右以及顶部过渡区。

在每一区域：

> 按 anchor angle / genomic order 排序。

然后沿外围逐行分配位置。

---

### 第五步：尽量保持原始顺序

这是避免连接线交叉的核心。

例如两个 anchor：

```text
A
B
```

在圆周上 A 在 B 前面，那么外面也尽量保持：

```text
A label
B label
```

而不是为了局部距离更短把两者交换。

---

### 第六步：统一避让

最小间距应同时作用于：

* restriction ↔ restriction
* feature ↔ feature
* restriction ↔ feature
* primer ↔ 任何其他 label

这样才能真正复制 pET-28a / pTRIPZ 的外观。

---

# 二十五、不同 annotation 保持各自样式，但共享空间

共享位置计算 ≠ 所有标签画成一样。

### Restriction site

裸文字：

```text
EcoRI (952)
```

灰 leader line。

### Feature callout

浅色圆角框：

```text
[ tet operator ]
```

浅色/灰 leader line。

### Primer

紫色文字：

```text
VR primer  (155 .. 174)
```

以及紫色 leader / binding marker。

它们只是：

> **共享 collision space 和排序规则。**

而不是共享同一种外观。

---

# 二十六、不要依赖“两个 repel 各退一点”解决问题

两个彼此独立的 repel engine 即使反复迭代，也很难得到 SnapGene 那种整齐结果。

原因是 SnapGene 的效果不是纯粹：

> “哪个碰到了就推开哪个”。

它还有：

* 圆周原始顺序；
* 左右区域；
* label 行；
* fan-out；
* leader crossing penalty；
* canvas expansion；
* map radius adjustment；

这些全局约束。

因此外侧标签应该被看作：

> **一个整体排版问题。**

---

# 二十七、当外围实在太拥挤时的处理顺序

为了接近 SnapGene，应优先：

1. 增大外围 label radius；
2. 让 label 向上下扩展；
3. 调整 map center；
4. 缩小 sequence circle；
5. 必要时扩大输出 canvas；
6. 才考虑隐藏低优先级内容。

而不是：

> 先缩字体。

SnapGene 的参考图中文字尺寸在非常拥挤时仍然保持可读，pETDuet-1 就是明显例子——它宁可把质粒圆画得很小，也没有把 enzyme labels 缩成难以阅读。

---

# 二十八、Primer 是你目前总结里缺失的一套视觉对象

这 13 个 `.dna` 还揭示了一个重要问题。

## pSB1C3

真正包含两个 Primer：

* VR primer
* VF2 primer

## pETDuet-1

真正包含四个 Primer：

* pET Upstream Primer
* T7 Terminator Primer
* DuetUP2 Primer
* DuetDOWN1 Primer

这些不是 `Feature type = primer_bind`。

例如 pBluescript 中：

`M13 fwd`

是：

**primer_bind feature**

颜色：

`#A020F0`

而 pSB1C3 的：

`VR primer`

则是 SnapGene 独立 Primer annotation。

图上真正 Primer 的视觉语言是：

* 紫色；
* binding location 在 DNA 附近用紫色短线/短括号表示；
* 标签可以移到圆外；
* leader line 也是紫色；
* 标签包含 range：

```text
VR primer  (155 .. 174)
```

所以如果 ggchord 将来支持 primer，应把它视为：

> **第三套 annotation object**

而不是偷偷转换成普通 feature。

---

# 二十九、Feature 与 restriction labels 的字体层次

整体应保持比较克制的字体体系。

### 中央质粒名称

最大、粗体。

### sequence length

略小、regular。

### feature labels

regular，通常不粗体。

### restriction enzyme name

unique → bold
nonunique → regular

### restriction coordinate

regular。

### primer

紫色 regular。

不要依赖很多不同字体大小制造层次，SnapGene 主要依靠：

* weight；
* color；
* position；
* feature geometry；

来建立信息等级。

---

# 三十、Feature label 与 arc 的关系

这也是当前 ggchord 与 SnapGene 很容易拉开差距的一点。

理想效果应是：

```text
       C a s 9
     ╭─────────→
```

而不是：

```text
       Cas9
     ╭─────────→
```

然后只把整个字符串旋转一个角度。

对于足够长的 curved feature：

> 各字符的局部 baseline 应逐渐跟随圆弧。

因此会产生 SnapGene 那种“文字天然长在箭头上”的感觉。

---

# 三十一、Feature 内外文字的位置不应机械绑定 feature midpoint

Feature 的几何 midpoint 可以作为首选 anchor，但不是最终 label 必须所在的位置。

需要考虑：

* feature 长度；
* 文本宽度；
* feature 前后是否存在其它 feature；
* 内层轨道；
* DNA tick；
* 内部其他 label；
* 圆外 restriction labels。

因此有时 label 会：

* 稍微偏离 feature 几何中点；
* 位于 feature 外侧邻近区域；
* 最终移动到圆外。

这种“先尝试贴附，再允许脱离”的行为比强制 midpoint 更接近 SnapGene。

---

# 三十二、非常短的 feature 不能因为 bp 少而消失

pET-28a、PX458、pEarleyGate 都能说明这一点。

例如：

* ATG 3 bp
* 6xHis 18 bp
* 小 peptide
* operator
* tag

即使实际 bp 极短，SnapGene 仍给它一个具有视觉存在感的 glyph。

因此应有：

> **minimum visual angular width / minimum visible glyph size**

但真实 genomic anchor 仍必须准确。

它解决的是视觉可见性，不是改变 biological coordinate。

---

# 三十三、配色的总体哲学

SnapGene 的成熟感并不来自“使用很多随机颜色”。

实际上它非常有规律：

### 高饱和主视觉

* yellow → origin
* burgundy → 重要 CDS
* pale green → selection/resistance
* teal → binding/recombination
* purple → primer
* bright green/red → fluorescent reporter

### 低饱和结构辅助

* white → promoter / terminator / enhancer
* gray → polyA / RBS / IRES/bom
* bisque → LTR/repeat/viral structural feature

因此：

> 调控骨架保持中性；功能主体获得颜色。

这点比简单复刻几个 HEX 更重要。

---

# 三十四、这批数据中还值得纳入 ggchord 的视觉规则

## 34.1 Origin seam 不能被当成真正断点

12 点只是坐标起点，不是物理断裂。

任何跨 origin 的：

* feature；
* label；
* promoter；
* MCS；

都要保持连续。

---

## 34.2 Feature outline 是重要组成部分

浅色 feature，尤其：

* `#FFFFFF`
* `#CCFFCC`
* `#FFFF00`

都需要明显但不过重的深灰细轮廓。

否则在 `#FFFCF5` 背景上白色 promoter 会消失。

---

## 34.3 空白也是 SnapGene 风格的一部分

不要追求：

> feature 填满整个圆内部。

pBR322、PX458 中都有大量空白。

Feature tracks 只在发生重叠时增加，不应该为了“对齐”而人为创建很多同心环。

---

## 34.4 同一类别重复出现必须保持一致视觉语义

例如 pTRE 的 7 个 tet operator：

* 同色；
* 同厚度；
* 同类 glyph；
* 同 label box 样式。

pETDuet 的两个 T7 cassette 同理。

这种重复非常适合检查 ggchord 是否出现“不稳定样式”。

---

# 三十五、数据库差异的验收原则

这点建议明确写进 Codex 任务。

**不要把“与 SnapGene 图上的 feature/site 数量逐项一致”作为验收条件。**

因为 ggchord 的：

* feature database；
* restriction-enzyme database；
* enzyme set；
* methylation metadata；
* common-feature detection；

都可能和 SnapGene 官方数据库不同。

允许出现：

* ggchord 多几个 feature；
* 少几个 feature；
* 多几个 restriction site；
* 少几个 restriction site；
* annotation name 略有不同。

这些都不是核心问题。

真正应该比较的是：

> 给定 ggchord 自己实际检测出来的 annotation 集合以后，其视觉呈现是否遵循 SnapGene 的规则。

重点比较：

**feature geometry
semantic color
overlap layering
arc labels
inside/outside decision
external callout style
restriction-label grouping
site fan-out
global collision avoidance
leader routing
map scaling
circular seam handling**

而不是 annotation 数目。

---

# 三十六、最终可作为视觉验收的核心规则

我建议 Codex 最终至少用这 13 个参考图检查以下事项：

1. **pBR322**：简单图不能因为算法复杂而变得难看。
2. **pUC19**：MCS restriction labels 必须形成整齐 fan-out。
3. **pBluescript**：长 feature 外、短 feature 内的 layering 必须正确。
4. **pSB1C3**：primer labels 与 BioBrick/terminator/enzymes 可以共存。
5. **pET-28a**：大量极短 feature 的外置 callout 必须和 restriction labels 统一避让。
6. **pETDuet-1**：极端标签密度下宁可缩小圆，也不能重叠；两个表达 cassette 保持一致。
7. **pcDNA3.1**：中等复杂度仍保持大面积空白和清晰层次。
8. **pTRE-Tight-BI**：双向 promoter、两个 MCS、7 个 operator 的结构关系清楚。
9. **PX458**：Cas9 数 kb 与数 bp feature 能在同一图中自然共存。
10. **pDONR221**：attP1/attP2/ccdB/CmR 的 Gateway 结构有明确语义色。
11. **pCAMBIA1300**：跨 origin feature 连续，T-DNA/bacterial backbone 清楚。
12. **pEarleyGate 201**：丰富的 feature 类型和颜色仍保持统一视觉语言。
13. **pTRIPZ**：作为最终 stress test，feature、site、外置 label、重复 operator、viral 元件均不发生无规则碰撞。

---

## 关于你第 6 点，我会把它定为这次修改的一个核心原则

**不要取消 `geom_feature_label_repel()` 和 `geom_restriction_site()` 两个公开入口；真正需要改变的是它们最终的“空间所有权”。**

现在如果是：

```text
feature layer
    ↓ 自己 repel

restriction layer
    ↓ 自己 repel

最后叠加
```

就天然存在互撞风险。

目标应该变成：

```text
features
   ├─ 能内部显示 → 内部
   └─ 不能内部显示 → external candidate
                              \
restriction sites ─────────────→ shared external layout
                              /
primers ─────────────────────/
```

最后才根据 annotation 类型分别应用：

* feature rounded callout；
* restriction plain enzyme text；
* primer purple text；

但它们的位置是**共同求解**出来的。

我认为这是你现在要从 SnapGene 学习的一个非常核心的架构性视觉规律：**Feature 与 restriction site 在数据语义上是两个系统，但在“圆外标签排版”上必须成为一个系统。**

尤其看 **pET-28a、pTRIPZ、pEarleyGate 201**，这一点几乎是直接写在图上的。只要 ggchord 继续让两个图层互不知情地独立 repel，即使单独看两个 geom 都很好，最终也很难达到 SnapGene 那种效果。

[1]: https://support.snapgene.com/hc/en-us/articles/10387315758100-SnapGene-Server-Request-API?utm_source=chatgpt.com "SnapGene Server Request API – SnapGene Support"
[2]: https://support.snapgene.com/hc/en-us/articles/10383910479124-Prioritize-Display-of-Features-in-Map-View?utm_source=chatgpt.com "Prioritize Display of Features in Map View – SnapGene Support"
[3]: https://support.snapgene.com/hc/en-us/articles/10383722725524-Display-Feature-Labels-Below-or-Inside-a-Map?utm_source=chatgpt.com "Display Feature Labels Below or Inside a Map – SnapGene Support"
[4]: https://support.snapgene.com/hc/en-us/articles/10383729806356-Don-t-Highlight-Unique-Cutters-in-Bold?utm_source=chatgpt.com "Don't Highlight Unique Cutters in Bold – SnapGene Support"
[5]: https://www.snapgene.com/updates/snapgene-5-1-3-release-notes?utm_source=chatgpt.com "SnapGene Version 5.1.3"
[6]: https://www.snapgene.com/updates/snapgene-version-8-2-0?utm_source=chatgpt.com "SnapGene Version 8.2.0"
