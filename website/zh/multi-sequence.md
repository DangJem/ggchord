## 同一数据，四种视图 {#task}

完整序列长度、真实 CDS 和已观测 HSP 可以共用四种关系坐标。[严格四视图对照](multi-sequence-views.html#input)保持相同输入、图层和映射，并显式显示全对关系。下图为独立设计的 Chord 旗舰；四张旗舰各自的显示选择和坐标适配见[旗舰说明](multi-sequence-views.html#flagship)。

{{manualfigure:flagchord}}

## 共用输入与图层 {#path}

[数据与图层](multi-sequence-start.html) · [四视图对照](multi-sequence-views.html) · [显示规则](multi-sequence-layout.html) · [连接](multi-sequence-links.html) · [标签](multi-sequence-labels.html)

## Chord：弧形关系 {#chord}
弧形骨架、半径、曲率和区间连接路由。[几何设置](chord-geometry.html) · [综合实例](chord.html)。

## Collinear：顺序与对齐 {#collinear}
横行或竖列排列、锚点、显示方向和局部窗口。[排列设置](collinear-layout.html) · [综合实例](collinear.html)。

## Radial：径向共线性 {#radial}
辐条序列、起始角度和角向排列。[径向设置与显示实例](radial-synteny.html)。

## Sankey：双列区间关系 {#sankey}
镜像 query/subject 列、列间距和重复显示实例。[双列设置](sankey.html)。

## 输入契约 {#inputs}
完整来源长度、标识和单位独立于显示几何。生物学链来自注释，切换坐标不会重写它。[数据契约](data.html) · [兼容范围](compatibility.html)。
