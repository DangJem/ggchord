已有注释只需长度与表格；本页为了获得真实示例注释才读取 DNA。

{{code:courseblockaj}}

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

零旋转的原点在 12 点钟，默认顺时针增加；正 rotation 逆时针转动整图。与上面 double 基线相比，下面分别只改变 direction 或 gap。rotation 对照见[圆图旋转](circular-layout.html#rotation)。

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

## 切换为单线骨架 {#style}

{{prepare:firstcircle,firstfeatures}}

{{code:frontcirclestyle}}

{{figure:frontcirclestyle}}

`seq_style = "single"` 把双线改为单线。中心名称和 bp 刻度仍可用于定位。

## 转动圆图 {#rotation}

{{code:firstfeatures,frontcirclerotation}}

{{figure:frontcirclerotation}}

默认 `rotation = 0` 把序列起点放在 12 点钟，坐标顺时针增加。这里设置 `rotation = 90`，起点移到 9 点钟；圆周刻度和复制起点一起旋转。原始注释区间仍是同一段 bp。下一节[整理自己的序列和注释](circular-data.html#tables)。
