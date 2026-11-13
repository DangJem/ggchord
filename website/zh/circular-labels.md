用相同六条真实 pUC19 注释比较固定与自动标签。位置、字号、换行与输出尺寸共同影响可读性；文字隐藏不代表注释消失。

## 独立准备本页 {#input}

需可选依赖 dplyr。原始名称保存在 finish_features，显示文字仅附加真实坐标，不修改注释区间或链。图层读取相同 finish_labels；下面所有图均为 9 × 8 英寸。

{{code:comparefinishcircle}}


## 固定标签：方向、侧向与旋转 {#fixed-labels}

基线 feature 方向跟随特征。后面依次只改为 horizontal、inside 或增加 30 度旋转。只改变一个设置时回到同一基线，不累计上张图的变化。默认 overlap 为 hide，因此拥挤文字可能隐藏。

{{code:comparefixedfeature}}
{{figure:comparefixedfeature}}
{{code:comparefixedhorizontal}}
{{figure:comparefixedhorizontal}}
{{code:comparefixedinside}}
{{figure:comparefixedinside}}
{{code:comparefixedrotate}}
{{figure:comparefixedrotate}}


## 固定入口的参数与冲突 {#fixed-parameters}

| 参数 | 默认 | 意义 |
| --- | --- | --- |
| label_orientation | feature | feature/horizontal/radial/tangent；水平与沿特征不同 |
| label_side | outside | outside/inside/auto |
| label_overlap | hide | hide/nudge/allow；allow 可能重叠 |
| label_rotation | NULL | 额外角度；默认解析为 0 |
| label_radial_offset | 0.04 | 局部法线方向的位移，不是 bp |
| label_circum_offset | NULL | 默认 0，沿周向位移 |
| label_circum_limit | NULL | 默认 TRUE，限制周向位移 |
| label_wrap | NULL | 可选字符换行宽度 |

这些固定控制不能传给自动 repel 入口。标签 Position 应与特征一致；子集/重排时优先使用唯一 feature_id，不能让同一 ID 指向不同区间。共线性中的 radial/tangent 对应竖排/横排，不沿圆弧排字。



## 同一文字比较四种自动布局 {#automatic-labels}

feature 默认先尝试内部，再尝试相邻与外部；radial、auto、callout 提供不同摆放策略。布局只改变显示，不能新增或删除源注释。四图保留相同输入、请求字号与尺寸。

{{code:comparefeaturefit}}
{{figure:comparefeaturefit}}
{{code:compareradialfit}}
{{figure:compareradialfit}}
{{code:compareautofit}}
{{figure:compareautofit}}
{{code:comparecalloutfit}}
{{figure:comparecalloutfit}}


## 换行、缩小、隐藏和引线 {#fitting}

从 feature 基线出发，先只将 label_fit 改为 shrink，再回到基线只指定 label_wrap=14。shrink 显式允许缩小，最低字号仍为默认 1.5；不能无限缩小来保证每个名称出现。wrap 的字符上限不是固定毫米宽度。

{{code:compareshrinkfit}}
{{figure:compareshrinkfit}}
{{code:comparewrapfit}}
{{figure:comparewrapfit}}

以 14 字符换行图为基线，只改 label_fit=ellipsis；随后回到 feature 默认图，只禁用 external。省略或隐藏只改变文字显示，不改源名称或特征。

{{code:compareellipsisfit}}
{{figure:compareellipsisfit}}
{{code:compareexternaloff}}
{{figure:compareexternaloff}}


## 自动入口的重要默认值 {#automatic-parameters}

| 参数 | 默认 | 选择与限制 |
| --- | --- | --- |
| label_layout | feature | feature/radial/auto/callout；gene 布局集合另见通用标签页 |
| label_fit | wrap | wrap/none/ellipsis/auto/shrink；省略文字不改源表 |
| label_max_lines | 2 | 最大行数；与 wrap 一起检查 |
| label_min_size | 1.5 | shrink 的字号下限 |
| label_side | outside | outside/inside/auto |
| external | TRUE | FALSE 时不使用外部回退；共线性不支持该约束 |
| max_overlaps | Inf | 默认不保证完全无碰撞 |
| label_segment_overlap | fade | 冲突引线的显示策略，默认淡化 |
| label_segment_overlap_alpha | 0.18 | 淡化后的透明度 |
| label_segment_linetype | auto | 自动选择引线线型 |
| label_segment_layer | below | below 置于数据层下方；above 恢复标签层顺序，文字位置不变 |

诊断表在 export_ggchord_layout 的标签记录中。先读隐藏/适配原因，再调整输出尺寸、标签子集、布局或换行。需要完整源名称时不要只保存截断后的图。




feature 布局使用直接引线，允许穿过特征，不再多折绕路。默认 below 时，正常特征
填充覆盖下方引线；透明填充仍保持透明。内部单行文字沿圆轨道排布，跨原点来源
只保留一条注释，闭合圆上的特征保持连续；有可见开口的坐标仍保留显示分片。相邻标签有明显角向位移时补短引线。
要恢复引线在标签层原有的绘制顺序，可使用
`geom_feature_label_repel(label_segment_layer = "above")`。

## 练习与自备注释 {#adapt}

只把 label_max_lines 改为 3，解释是否改善了拥挤；再固定输出尺寸，比较外部回退开关。固定标签与自动标签有各自接口，不能混用。自备数据用明确 label 文本或 anno，核对 ID、区间、链及唯一 feature_id。继续看[引物文字](primers.html#primer-labels)和[密集 pETDuet-1 案例](plasmids.html#complex)。


相邻注释优先比较实际径向距离与角向位移，内侧有空闲轨道时不强迫它们在同一
最近轨道上挤开。箭头中的虚线可能是来源分段或 cleavage 边界，应检查返回的
`segments`、`segment_name` 和 `cleavage_arrows`，不能把每条虚线当作独立特征。
顶层特征标签不会自动为每个子段另加文字；需要子段名称时保留并查看源分段表。
