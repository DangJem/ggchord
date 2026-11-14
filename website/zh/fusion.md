## 独立视图与来源区间连接 {#task}

Fusion 是组合层。综合示例在明确的面板矩形中放置 Chord 总览、Collinear 视图、Radial 视图和双序列局部细节。连接来自已观测 HSP 和注册来源区间；同一记录的重复显示不构成新增观测。

{{manualfigure:flagfusion}}

## 组合与布局 {#path}

[面板与来源连接](fusion-start.html) · [位置、尺寸与图例](fusion-layout.html) · [跨面板路由](fusion-routing.html)

## 局部放大 {#magnify}

以明确选区连接全局图和独立构造的局部图；注册的真实对齐字符可进一步提供残基或碱基细节。见[局部放大与字符细节](fusion-details.html)。

## 嵌套与环状选区 {#nested}

嵌套 Fusion 保留叶面板路径；跨原点选区需要显式环绕声明和匹配的完整来源长度。见[嵌套组合与跨原点放大](fusion-advanced.html)。

## 支持的面板 {#inputs}

按各自契约接受注册基因组视图、嵌套组合和注册 Cartesian 字符细节。Dotplot、分面和未注册外部图不做自动猜测或适配。
## 旗舰图的展示选择 {#flagship-selection}

旗舰图的三个概览面板各保留32条原始HSP和266条CDS；面板之间仅选3条真实观测展示对应关系。局部图选择长度至少2000 bp的观测中最短的一条，并保留两侧各800 bp上下文（在序列边界截断），展示原有CDS及名称。淡色跨面板带标记对应区间，面板内部的ribbon颜色表示Identity；不要把展示选择或重复视图当成新增观测。
