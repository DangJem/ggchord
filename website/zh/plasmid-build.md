这一页只关注每次新增的视觉信息。准备时计算输入，绘图仍由独立图层完成。

## 本页准备 {#prepare}
{{prepare:featurepalette,docstyle,plasmidinput}}
{{code:learnplasmid}}

## 第一步：双线骨架与中心名称 {#backbone}
{{figure:plasmidbackbone}}
先检查长度和名称。`coord_circular()` 仅接受一条序列；双线样式由geom控制。

## 第二步：加入特征形状 {#features}
{{figure:plasmidfeatures}}
`position_feature_stack()` 为重叠注释分层，`geom_feature_plasmid()` 按实际类别选择形状和配色。短特征可能使用紧凑形状，不等于缺失。

## 第三步：为特征适配标签 {#labels}
{{figure:plasmidtext}}
箭头层与标签层使用同一Position。允许文字在明确字号下限内缩小，必要时外置；尺寸变化后重新检查。

## 第四步：引物与选定酶切位点 {#details}
{{code:plasmid}}
{{figure:plasmid}}
只展示选定的三个酶切位点，完整搜索表仍保留。完整的识别方法、来源与假设见[3.3](recognition.html)。然后使用[13个参考案例](plasmid-catalog.html)比较注释密度，最后完成[pETDuet-1综合图](plasmids.html#complex)。
