<div class="task-intro"><p>从已知 DNA 搜索引物、限制酶和常见注释，再决定如何显示；计算命中与实验验证是不同事情。</p></div>

## 本页准备 {#prepare}
本例使用随包 pUC19。

## 搜索与显示分开 {#search}
直接使用随包 pUC19 DNA，查找目录引物的精确命中与选定酶切位点：

本例整理序列表时使用 dplyr。先运行 `library(dplyr)`。

{{code:recognition}}
{{figure:recognition}}

圆形搜索也检查跨原点的候选位置。线性序列应设置 `circular = FALSE`。搜索命中是计算匹配，
不等于实验验证。引物 `"three_prime"` 模式寻找精确后缀，没有错配模型；
Wallace 熔解温度估计是基于组成的启发式值。

## 解释酶切结果 {#restriction}
`find_restriction_sites()` 返回未筛选匹配和切割坐标，`filter_restriction_sites()`
选择显示集合。酶位点计数针对完整输入序列，筛选后保持不变。
`methylation = "auto"` 仅为精确识别的内置参考使用已审计背景，不能据此声称未知
样品不受阻断。本例明确选择 `"none"`，这是搜索条件，不是对实际样品甲基化状态的测量。
`geom_restriction_site(leader = "radial")` 是默认连线方式。

## 通用注释识别 {#features}
`find_common_features()` 支持精确参考的整理注释、DNA 匹配及逐记录检测模式。
蛋白搜索采用标准遗传密码六框翻译，不处理剪接、替代起始或跨原点蛋白匹配。
来源中的遗传密码元数据不选择另一套翻译表；非标准遗传密码应使用 DNA 匹配或外部
支持相应遗传密码的分析。confidence 是启发式标签，不是概率。
已知蛋白/基因组对应关系应使用[显式 CDS 投影](projection.html)。

## 跟随完整参考注释流程 {#reference}
这个简短流程和完整圆图使用同一条有来源的随包质粒。[完整 pUC19 流程](plasmids.html)将目录注释、
精确参考匹配、酶位点显示筛选、引物标签及注释堆叠组合起来。
只修改序列名称不能让它成为已知参考，有来源的识别依据是序列内容。
应读取结果方法及甲基化来源，不能把参考结论赋予未知样本。

## 检查要点 {#check-result}
查看命中坐标、方法与来源。酶位点筛选是显示选择；甲基化搜索条件不是对未知样本状态的测量。

<div class="next-step">继续：<a href="primers.html">识别引物、弧线与文字 →</a></div>

识别结果确认后，分别复现[引物弧线与文字](primers.html)和[酶切标签与引线](restriction-sites.html)。圆图骨架、中心文字和分层见[图层课](circular-layout.html#input)，固定与自动文字见[圆图标签](circular-labels.html)。

## 目录不是搜索结果 {#catalog}

primer_catalog 返回寡核苷酸目录；aliases=TRUE 返回来源别名。目录方向是描述提示，实际结合方向由 find_primer_bindings 在目标 DNA 上计算。同序列多名称不应被当成多个不同寡核苷酸。按需查看目录后，去[引物课](primers.html)搜索和显示。

{{code:comparecatalog}}
{{console:comparecatalog}}
