<div class="task-intro"><p>从已知 DNA 搜索引物、限制酶和常见注释，再决定如何显示；计算命中与实验验证是不同事情。</p></div>

## 本页准备 {#prepare}
在新的 R 会话中先运行下方代码，再按页面顺序继续。配色和排版是可选的显示设置，输入只选择需要的列。

```r
library(ggchord)
library(dplyr)
```

## 搜索与显示分开 {#search}
直接使用随包 pUC19 DNA，查找目录引物的精确命中与选定酶切位点：

{{code:recognition}}

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

## 完成后检查 {#check-result}
查看命中坐标、方法与来源。酶位点筛选是显示选择；甲基化搜索条件不是对未知样本状态的测量。

<div class="next-step">继续：<a href="plasmids.html">绘制完整质粒图 →</a></div>
