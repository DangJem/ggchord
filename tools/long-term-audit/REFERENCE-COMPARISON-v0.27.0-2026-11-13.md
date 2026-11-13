# 13张质粒参考图逐项比较（2026-11-13）

已逐张查看 examples/plasmid 的13张原PNG，独立读取相应.dna的Feature、Segment、
Primer和酶profile，与当前完整find_common_features/find_restriction_sites结果核对。
不是OCR推测，也不以缺少导出设置为开发阻塞。原始PNG不改。

## 可以确认的来源差异

- 201条顶层Feature全部识别且一一匹配：名字、类型、start/end/strand无差异。
  完整源表和识别表见reference-source-results.json。不能把这项数据核对等同于
  全部文字均显示、显示形状正确或PNG与.dna导出设置完全一致。
- 旧图库/用户截图每张只选择最靠近三个等分目标的3条识别记录，不是只有3处位点。
  完整REBASE包含大量非商业/同识别谱酶，不应把所有原始记录都塞进参考图。
  本轮比较明确使用保存profile及commercial_nonredundant父集合；同位点酶名可合并。
- pSB1C3的2条、pETDuet-1的5条独立primer不属于顶层Feature，旧图库没有加primer层。
  参考图中的M13等有时本身是primer_bind Feature，不能把两种来源简单合并计数。
- 10个来源含17个具名子段；旧数据生成器把名字和translated都写为缺失/假。
  已修生成规则并实际重建，恢复17个名字和64个translated标志。逐对象比较证实
  REBASE/primer数据库及common-feature其余字段不变。没有改原始DNA或推断新名字。
  原子段表单独导出到named-segments.csv。顶层标签默认仍只有一条，不自动给每个
  子段生成独立注释；“元数据修复”与“子段文字显示”分别报告。

箭头内部虚线可能表示分段边界、cleavage标记，或源记录显式线型，不是统一的
“特殊区域”。例如pTRE AmpR：1541–2332主段、2333–2401具名signal sequence，
并记录cleavageArrows=2332。参考PNG也有这段边界及截短的si...字样。lac promoter
的-35/-10是具名子段。短片段是否另加标签需要以真实来源和空间为依据，不能从
虚线自行创造新的顶层feature。

## 逐图观察与记录

下表数量是.dna顶层Feature/原独立primer，以及本次明确策略选中的酶–位点记录数，
不是人工数PNG标签，亦不是独立切割位置数。位置数与完整记录另见CSV。

| 参考图 | Feature/独立primer | 当前profile记录数 | PNG与源/绘图需要区分的点 |
| --- | --- | --- | --- |
| pUC19 | 9/0 | 39 | AmpR信号子段、lac promoter的-10/-35；MCS/M13在内层；原图右上MCS位点扇束，旧ggchord抽样遗漏了展示密度 |
| pBR322 | 6/0 | 47 | 顶层6特征一致；AmpR分段虚线；原图部分位点有灰色星号，依赖甲基化语义；同位点等效酶名称与REBASE代表不同 |
| pSB1C3 | 7/2 | 38 | VR/VF2是独立引物，旧图缺少两层；顶部终止子/BioBrick相邻拥挤，原图也将一条外引，不应误认为未识别 |
| pBluescript II SK(+) | 14/0 | 48 | 顶层含primer_bind而无独立Primer；原图部分SK primer未显字，ggchord全源显示可多出字；AmpR及lac promoter具名子段 |
| pET-28a(+) | 16/0 | 51 | 原图MCS/ATG/RBS/tag等多条外引；ggchord优先内轨，短特征空间不足需算法处理；酶等效命名不同 |
| pETDuet-1 | 20/5 | 保存profile为0 | 原PNG有大量位点和独立primer，.dna却保存None；比较图明确另选unique_6plus，不伪称恢复原图设置；跨原点T7需连续 |
| pcDNA3.1(+) | 13/0 | 51 | CMV enhancer/promoter邻近，AmpR有signal sequence；原图箭头连接/分段与位置层级不同，不能用顶层一致证明外观一致 |
| pDONR221 | 11/0 | 47 | M13等为Feature，attP/终止子字多；原图外引M13 fwd，ggchord内轨；切点等效名称、甲基化星号需分开 |
| pCAMBIA1300 | 17/0 | 42 | lacZα及MCS跨原点；源lac promoter含-35/-10；左侧LB标签应沿内轨；原图多数位点在顶部MCS区域 |
| pTRE-Tight-BI | 15/0 | 47 | 跨原点双向promoter的假接缝已修；7个operator均有来源；左上长字挤同轨需改候选评分；AmpR信号段有名但旧导入丢失 |
| pEarleyGate 201 | 21/0 | 39 | 原PNG右侧多个短特征外引；lac UV5 promoter有-35/-10；原图LB外引，不能把标签层差异当源缺失 |
| pTRIPZ | 33/0 | 38 | 保存Unique Cutters+BamHI自定义集合；原图多条operator/TurboRFP外引、IRES等邻接密集；AmpR含信号段，不以这张图进行生物工程解释 |
| PX458 | 19/0 | 4 | 自定义BbsI+EcoRI恰4处，原图少位点是保存选择；多色Cas9/EGFP/NLS等保持独立源，EGFP 1a与AmpR信号子段需保留 |

当前REBASE代表酶不一定等于旧PNG名称，例如同识别谱的代表替换、供应商集合变化、
同位点多个不同motif、nickase与双链切割语义都会改变显示列表。没有从PNG推断未
提供的全部酶库版本/供应商选择；saved profile不是完整导出记录。甲基化假设和
被阻断位点仍保留自己的字段，不把灰色星号直接翻译成“没有识别”。

## 比较图和布局范围

`render-reference-catalogue.R` 是本轮独立比较脚本，保留全部顶层feature、保存的
独立primer及上述酶筛选，10×8in/110dpi，不覆盖原图库/冻结网站。相邻标签评分的
第一候选lane图保留原3位点同输入尺寸控制。源表核对不依赖后续PNG质量通过；
每张图的可见/隐藏、物理收敛、警告和冲突都保留到独立results，未通过不清零。

本次没有进行发布级测试、全设备/字体或SnapGene同尺寸等价认证，也没有关闭
长期任务书的全部布局质量目标。


## 候选迭代与运行边界

第一轮全量位点候选在完成4张后，pET-28a运行持续占用CPU，定位为旧限制位点
多折graph避障后主动发SIGINT停止。退出1和部分结果保存为rejected-partial，
不是完整13图成功。随后去除圆注释的多折graph后处理；保留外圈队列原本的直线/
单shoulder路径，最终接触仍计数，不以绕远掩盖冲突。

内轨评分第一候选过度惩罚小角向移动，在pCAMBIA让M13等名字形成较深径向列。
最终候选允许约四分之一字高的局部调整，再比较角向位移与真实径向距离；
ordered packing不再具有无条件的负分优先权。失败候选lanes文件保留，不算
最终视觉通过；balanced文件为后续最终两例稀疏控制，不用它替代密集图指标。


## 最终本轮实际结果

13张最终PNG已完整退出0，源快照校验通过，均人工查看。全部757个源标签均可见，
文字重叠/画布外标签为0；12张物理收敛，pET-28a未收敛，其67标签仍可见但测量
属于估计。这不是全部布局质量通过。以下保留未消除的中心线路径接触：

| 案例 | 可见/源标签 | 收敛 | 引线交叉 | 穿字框 |
| --- | --- | --- | --- | --- |
| pUC19 | 48/48 | 是 | 0 | 2 |
| pBR322 | 53/53 | 是 | 0 | 0 |
| pSB1C3 | 47/47 | 是 | 0 | 5 |
| pBluescript_II_SK_plus | 62/62 | 是 | 1 | 3 |
| pET_28a_plus | 67/67 | 否 | 9 | 17 |
| pETDuet_1 | 83/83 | 是 | 1 | 12 |
| pcDNA3_1_plus | 64/64 | 是 | 0 | 0 |
| pDONR221 | 58/58 | 是 | 0 | 1 |
| pCAMBIA1300 | 59/59 | 是 | 2 | 8 |
| pTRE_Tight_BI | 62/62 | 是 | 2 | 9 |
| pEarleyGate_201 | 60/60 | 是 | 1 | 9 |
| pTRIPZ | 71/71 | 是 | 1 | 3 |
| pSpCas9_BB_2A_GFP_PX458 | 23/23 | 是 | 1 | 10 |

同位点合并文字会保留多个来源表示，源标签总数不是文字块数。fade下中心线
穿框不等于最终可见墨迹穿框；本轮没有逐墨迹遮挡认证。pET-28a密集top-fan的
计算及非收敛也未修复，取消多折graph不是所有性能问题的完整根因结论。
最终比较图保留约9分钟运行中的真实失败，不以更大画布偷偷替换原10×8设置。

两张balanced稀疏控制均收敛、全部38标签可见，无自身弧字溢出/偏移无引线/
特征引线超一拐点/内部单行水平回退；AmpR promoter与SV40 signal分到不同内轨。
图内仍可见一些长径向引线和紧密小字，整体现阶段阅读质量仍待改善。

随包13条DNA字符串与reference逐条完全相等，独立记录见bundled-reference-equality。
最终render/source/balanced三个报告逐R SHA核对一致；sysdata SHA另由segment-
metadata-results保存，未把旧13目录/整站结果升级为当前布局认证。原201顶层
特征一致、17具名子段和7独立primer的事实与这些密集图质量失败分别保留。
