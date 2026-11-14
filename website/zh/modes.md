## 四大模块及其职责 {#modules}

| 模块 | 主要用途 | 坐标系或组合层 |
| --- | --- | --- |
| [模块一 · 多序列关系](multi-sequence.html) | 比较序列区间、注释和连接 | Chord、Collinear、Radial、Sankey |
| [模块二 · 单序列圆图](circular.html) | 在一条完整序列上组织注释 | Circular |
| [模块三 · HSP 点图](dotplot.html) | 以双轴显示 query/subject 的 HSP 端点 | Dotplot |
| [模块四 · 多视图组合](fusion.html) | 排列独立注册图并连接真实源区间 | Fusion 组合层 |

六种 Coord 负责单张图的全局坐标；Fusion 组合独立图形，不是第七种坐标系。Dotplot 使用专属线段图层，目前不作为 Fusion 的注册端点。

## 共用的绘图语法 {#prepare}

| 组成 | 职责 | 入口 |
| --- | --- | --- |
| 数据 | 来源标识、完整长度、区间与单位 | `ggchord()` |
| 几何图层 | 骨架、连接、特征、轨道和文字 | `geom_*()` |
| 映射 | 将数据字段关联到视觉属性 | `aes()` |
| 标度 | 规定颜色、取值范围和图例 | `scale_*()` |
| 局部位置 | 安排特征、标签和轨道 | `position_*()` |
| 坐标系 | 排列整张图的序列和区间 | `coord_*()` |
| 主题 | 控制非数据外观 | `theme_*()` |

## 选择关系视图 {#choose}

### Chord：弧形关系总览 {#chord}
沿弧形骨架组织多序列，便于观察整体连接。见[骨架、曲率与半径](chord-geometry.html)。

### Collinear：顺序与局部对齐 {#collinear}
将序列排成横行或竖列，便于比较顺序、锚点与局部窗口。见[对齐与排列](collinear-layout.html)。

### Radial：径向共线性 {#radial}
将序列作为辐条，保留各自坐标并在中心区域展示共线性。见[径向布局设置](radial-synteny.html)。

### Sankey：双列区间关系 {#sankey}
每条来源序列在 query 和 subject 两列各显示一次；连带表达观测区间，不代表守恒流量。见[双列布局设置](sankey.html)。

### Circular：单序列注释 {#circular}
一条序列构成圆环，分别叠加已有注释、引物和限制位点。见[单序列圆图](circular.html)。

### Dotplot 与 Fusion {#composition}
Dotplot 在两条序列轴上连接 HSP 四端点。Fusion 排列已经独立有效的图，并连接明确的来源区间。两者承担不同职责，不能互换图层契约。

## 来源身份与单位 {#terms}
生物学链、显示方向和比对方向分别记录。bp 与 aa 属于不同坐标域，跨域投影需要明确的 CDS 映射；完整序列长度不能从最后一个可见特征推断。

{{terms}}


<span id="check-result"></span><span id="coordinate-controls"></span>
