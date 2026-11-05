在同一条真实 pUC19 上调整骨架、中心文字和注释位置。先确认起点与方向，再添加信息层；需要的是单条完整 DNA 及有来源的注释。

## 本页输入与来源 {#input}

本页在新 R 会话独立运行，需可选依赖 dplyr。`finish_features` 保留识别结果的前六行，不改变区间与链。为比较标签适配，`finish_labels` 的显示文字在原始名称后附上真实坐标，原始表仍单独保存。

{{code:comparefinishcircle}}


## 骨架样式与双线间距 {#backbones}

`seq_style = "auto"` 在单圆图中选择双线；single 是一条骨架，double 是双线。下面保持 9 × 8 英寸与同一序列，每次只改变样式或双线间距。双线间距 `seq_backbone_gap` 默认 0.025，band 的宽度参数 `seq_backbone_width` 默认 0.035，二者不是线条的 `linewidth`。

{{code:comparebonedouble}}
{{figure:comparebonedouble}}
{{code:comparebonesingle}}
{{figure:comparebonesingle}}
{{code:comparebonegap}}
{{figure:comparebonegap}}


## band 色带骨架 {#band-rendering}

band 用填充多边形表示整条骨架，宽度由 seq_backbone_width 控制，默认 0.035。以下使用真实 band 骨架，与双线图相比只改变 seq_style。

{{code:compareboneband}}
{{figure:compareboneband}}

同一 band 基线，只把宽度加倍；其他设置不变。

{{code:comparebandwidth}}
{{figure:comparebandwidth}}

## 圆图坐标、旋转、开口与适配 {#coordinates}

开发版零旋转的原点在 12 点钟，默认顺时针增加；正 rotation 逆时针转动整图。与上面 double 基线相比，下面分别只改变 direction 或 gap。rotation 对照见[圆图旋转](plasmid-tuning.html#rotation)。

| 参数 | 默认 | 类型与边界 |
| --- | --- | --- |
| `rotation` | 0 | 有限数值，角度；90 是九点钟 |
| `direction` | clockwise | clockwise/counterclockwise；不是生物学链 |
| `gap` | 0 | [0,360) 度；不从序列中删除碱基 |
| `fit` | labels | labels/geometry/manual；是否将标签纳入适配 |
| `ratio` | 1 | 全局纵横比，不改变长度单位 |
| `clip` | off | 裁剪由坐标控制；打开裁剪可能截掉外侧文字 |

固定尺寸下，自动适配会改变整体缩放和可读刻度标签的子集。手动限制 `xlim/ylim` 时先检查文字和引线，不把裁剪误认为数据丢失。圆图用 gap/direction，不能传入 geom_seq 的 seq_gap/seq_orientation。

{{code:comparecirclecounter}}
{{figure:comparecirclecounter}}
{{code:comparecirclegap}}
{{figure:comparecirclegap}}


## 分别设置中心名称和长度 {#center}

中心文字与沿骨架的序列标签是不同图层。`name = TRUE`、`length = TRUE`，默认用换行符 separator 分隔；name 也可提供显示字符串。name_style/length_style 默认空列表，接受固定 colour、size、family、fontface、lineheight。下面先只隐藏 length，再回到基线只改变 name_style 的 size，长度字号不会跟着改变。

{{code:comparecenter}}
{{figure:comparecenter}}
{{code:comparecenternolength}}
{{figure:comparecenternolength}}
{{code:comparecenterstyle}}
{{figure:comparecenterstyle}}


## 通用形状与质粒预设 {#feature-shapes}

先用同六条真实注释只改变 arrow/block。block 不画方向箭头，但原始 strand 没有改变。第三幅展示 geom_feature_plasmid 的组合预设，默认 width=0.058、arrow_head_length=0.068、arrow_head_width=1.42，并使用质粒 Position；这会同时改变多个默认，不能当作只变一个参数。其他形状用途见[形状表](annotations.html#shapes)，箭头细节与短特征见[几何课](feature-geometry.html)。

{{code:comparefeaturearrow}}
{{figure:comparefeaturearrow}}
{{code:comparefeatureblock}}
{{figure:comparefeatureblock}}
{{code:comparefeaturepreset}}
{{figure:comparefeaturepreset}}
## 让重叠特征与文字共享位置 {#stacking}

下面用同一六条实际注释比较 identity 与共享的 feature stack。Position 同时传给几何和文字，不单独挪动文字假装注释分层。`position_feature_stack(spacing = 0.1, side = "strand", base_position = NULL)` 先应用可选基础 Position，再分配重叠轨道。`position_plasmid()` 是局部预设；`position_strand()` 按真实链应用带符号位移。检查方向时不要用内外侧代替生物学链。

{{code:comparefinishstack}}
{{code:comparepositionidentity}}
{{figure:comparepositionidentity}}
{{code:comparepositionstack}}
{{figure:comparepositionstack}}


## 练习、替换输入与下一步 {#adapt}

只把双线间距改为 0.04；再把中心文字的 length_style 改为不同字号。保留完整长度，不用注释最大 end 猜长度。用自己的 FASTA/GenBank/GFF3 时核对 ID 和坐标，见[格式导入](formats.html)。接下来学习[圆图标签](circular-labels.html)、[引物](primers.html)、[酶切](restriction-sites.html)和[轨道](tracks.html)。所有预览使用图中尺寸；用同样宽高的 ggsave 保存。
