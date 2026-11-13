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
| label_repeat | individual | bracket 为连续同名短特征添加共享的名称×数量与弧形括注 |
| label_segment_colour | auto | 圆内特征连圈外标签默认灰色；可指定颜色 |

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


密集内圈注释按来源角度顺序处理，候选位置同时检查文字与引线。引线可穿过
特征，但不为避让其他文字增加多次拐弯；优先另选合法文字位置。启用external
时，外层重复短注释可以成组使用外圈；较内层的特征先尝试内圈，空间不足时仍可回退。该策略是有界搜索，
仍应查看导出中的隐藏、收敛和冲突记录，不保证任意密度都能排入内圈。


## 连续重复特征的共享标注

`label_repeat = "bracket"` 是显式的显示样式，默认 `"individual"` 仍逐条标注。
同一特征轨道内，相邻、不重叠、名称、类型、方向和已解析填充一致的短特征，
可共享一个沿弧的“名称 × 数量”和弧形范围括注。中间插入其他注释时不跨过它
合并。各个箭头、坐标和源记录不合并；导出的可见性与适配记录将这些来源标为
`represented_by_shared_label`。此样式只用于单序列闭合圆的 `feature` 布局，
其他坐标保持逐条显示。弧形括注是范围标记，标签引线仍使用简单路径。

```r
features <- find_common_features(plasmid_example_pTRE_Tight_BI)
feature_position <- position_feature_stack(base_position = position_plasmid())
ggchord(plasmid_example_pTRE_Tight_BI[c("accver", "length")], feature_data = features) +
  geom_seq(seq_style = "double") +
  geom_feature_plasmid(position = feature_position) +
  geom_feature_label_repel(
    position = feature_position,
    label_repeat = "bracket",
    label_fit = "shrink",
    label_min_size = 2.5
  ) +
  coord_circular() +
  theme_ggchord_plasmid()
```

质粒主题的自动轴文字贴近序列骨架，显式轴偏移不变。圆形特征轨道按不能放入
箭头的文字需求留空，文字轨道按实际字号与物理间隙排布；显式 `spacing` 仍是
最小间距。这是有界布局，仍需按实际导出尺寸查看诊断，不能保证任意密度无冲突。


默认灰色虚线用于**从序列圆内的特征连到圈外标签**的引线。圈内相邻标签及限制
位点引线保留原有风格。若圈内名字紧贴自己的特征、角向未明显偏移且无竞争来源，
省略无助于归属的短线；偏移或邻近其他来源时仍保留引线。

```r
geom_feature_label_repel(
  label_segment_colour = "grey50",
  label_segment_linetype = "dashed"
)
```

显式参数调整该图层全部特征引线；自定义 `ggchord.feature.label.segment` 主题也
优先于自动配色。分段虚线按实际箭头轮廓裁切，保留原分段信息。圈内换行文字的
各行沿各自同心圆逐字排布。自动轴文字保留小的物理留白，显式偏移值仍优先。
