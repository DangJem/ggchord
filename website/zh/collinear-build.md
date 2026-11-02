先验证共同语法，再为直线视图调整排版。以下三个对象只在弦图教学对象上替换坐标系，不改表格、不重新计算比对。

## 本页准备 {#prepare}
先解压[教学输入](../assets/tutorial-inputs.zip)。
{{prepare:sequencepalette,docstyle,showcasedata}}
{{prepare:learnseq,learnlinks,learngenes,learnlabels}}

## 替换coord {#switch}
{{code:learnlinear}}

## 从骨架看比例 {#backbones}
{{figure:linearseq}}
每一行共享同一bp比例。序列长度不同，行长也不同；不是将每行任意缩放成一样长。

## 加入连接 {#links}
{{figure:linearlinks}}
相同的HSP保留在图层中，跨行连接也保留，因此可直接看到默认路由的局限。长连接更适合在下一页明确路由。

## 加入上下分链与标签 {#labels}
{{figure:linearlabels}}
箭头和标签来自同一个筛选表。正负链不改写；显示翻转也不会更改原始strand。继续：[2.3 调整跨行路由](collinear-routing.html)。
