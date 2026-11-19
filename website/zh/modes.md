<span id="input"></span><span id="read"></span><span id="save"></span><span id="source"></span><span id="stable"></span><span id="verified-versions"></span>
<a id="modules"></a>

## 想画哪一种图？ {#modules}

| 我想…… | 查看指南 | 直接试画 |
| --- | --- | --- |
| 比较多条序列的顺序、注释和区间连接 | [多序列关系](multi-sequence.md) | [从序列骨架开始](multi-sequence-start.md#draw) |
| 沿一条环形序列显示特征和轨道 | [单序列圆图](circular.md) | [从长度和已有注释开始](circular-start.md#input) |
| 在两条坐标轴上查看局部比对 | [序列比对点图](dotplot.md) | [画出 HSP 的端点和方向](dotplot-start.md#dotplot) |
| 组合总览、局部放大或不同类型的图 | [多视图组合](fusion.md) | [组合并连接两张图](fusion-start.md#panels) |

还可以在[案例图库](gallery.md)按图形和绘图任务查找可复现的例子。

<a id="choose"></a>

## 多序列关系选哪种布局？ {#choose}

<a id="chord"></a>

### Chord · 弧形关系总览 {#chord}

沿弧形骨架显示多条序列及其连接。查看[骨架、曲率与半径](chord-geometry.md)。

<a id="collinear"></a>

### Collinear · 顺序与局部对齐 {#collinear}

将序列排成横行或竖列，适合比较基因顺序、锚点与局部区域。查看[对齐与排列](collinear-layout.md)。

<a id="radial"></a>

### Radial · 径向共线性 {#radial}

把序列排成辐条，保留各自坐标，在中心区域显示共线关系。查看[径向布局](radial-synteny.md)。

<a id="sankey"></a>

### Sankey · 双列区间关系 {#sankey}

把序列分到 query 和 subject 两列，连带表示实际观测的区间对应。查看[双列布局](sankey.md)。

<a id="circular"></a>

### Circular · 单序列注释 {#circular}

沿一个圆环排列一条序列的特征、引物、限制位点或数值轨道。查看[单序列圆图](circular.md)。

<a id="composition"></a>

### Dotplot 与 Fusion {#composition}

[Dotplot](dotplot.md)在两条序列轴上显示 HSP 的成对端点。[Fusion](fusion.md)排列独立图形，可通过它们所表示的源区间添加跨面板连接。

<a id="shared-techniques"></a>

## 图画出来后，怎样调整？ {#shared-techniques}

| 目标 | 操作 |
| --- | --- |
| 调整箭头形状与放置方式 | [特征形状](feature-geometry.md)和[重叠注释](annotations.md) |
| 处理拥挤标签 | [文字适配与避让](labels.md) |
| 添加信号、覆盖度或其它数值轨道 | [轨道](tracks.md)和[轨道统计](tracks-statistics.md) |
| 设置配色、图例或主题 | [图层样式](styling.md)和[标度](scale-families.md) |
| 按序列或链分别设置几何参数 | [参数写法](parameters.md) |
| 按指定尺寸保存，或检查布局问题 | [输出](output.md)和[诊断](diagnostics.md) |

<a id="prepare"></a>

## 同一套绘图语法 {#prepare}

| 组成 | 职责 | 入口 |
| --- | --- | --- |
| 参考数据 | 提供序列 ID、长度和单位 | `ggchord()` |
| 几何图层 | 绘制骨架、连接、特征、轨道和文字 | `geom_*()` |
| 映射 | 把数据字段关联到视觉属性 | `aes()` |
| 标度 | 规定取值如何显示及图例如何呈现 | `scale_*()` |
| 局部位置 | 安排特征、标签和轨道 | `position_*()` |
| 坐标系 | 决定一张图的整体布局 | `coord_*()` |
| 主题 | 设置字体、间距和其它非数据外观 | `theme_*()` |

<a id="terms"></a>

## 来源身份与单位 {#terms}

序列的来源标识、坐标单位和完整长度来自输入。同一来源可以有多个显示实例，显示方向也可独立调整。核酸的 bp 与蛋白的 aa 属于不同坐标域，跨域请使用[明确的坐标投影](projection.md)。

{{terms}}

<span id="check-result"></span><span id="coordinate-controls"></span>

<a id="getting-started"></a>

## 开始使用 {#getting-started}

[安装并画出第一张图](installation.md#first-plot)。有自己的表格时，可直接查看[最小字段与坐标约定](data.md#tables)。
