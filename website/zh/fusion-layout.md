{{code:courseblockbm}}

{{prepare:fusiondata,fusionpanels,fusionlinear,fusionview}}

## 样式与物理布局 {#styles}

跨面板层接受`aes(fill=..., colour=..., alpha=..., linewidth=..., linetype=...)`及标准ggplot2标度。子图的角色标度只在各自面板生效；明确指定的固定样式优先于对应映射。默认Identity范围0–100，按0/50/80/90/95/100分段并保持数值比例。端点渐变与subject尖头须显式开启，复用单图Ribbon绘制器。尖头是装饰性方向标记，仍导出完整来源前缘。

{{code:fusionstyle}}

默认`layout="fixed"`分配面板区域，并按输出尺寸重新测量标签和图例；`layout="adaptive"`还按物理宽高比调整中央区域，并对外层原生Guide联合预算。用户明确指定的字号保持固定。这不是通用字号缩放或全局碰撞求解器。
每张子图保留自己的图例，不自动合并含义不同的标度；需要简化时，主动设置`show.legend=FALSE`或局部`theme()`。组合层的`guides()`配置跨面板Guide；`theme()`还配置后文的原生外层图形装饰。必须保留的图例空间会实际测量；画布过小则给出错误诊断。
## 八方位与独立物理区域 {#advanced-placement}

`fusion_panel()` 在四个基本方向之外支持 `topright`、`topleft`、`bottomright`、`bottomleft`。默认角落布局为八个方向预留空间，每个侧边名称只使用一次。其他排布可用 `angle`/`distance` 或显式 `box`，显式区域/角度没有固定面板数上限，但每个区域仍须可容纳。位置参数只改变容器，不旋转子图Coord或文字。角度在矩形画布上也表示真实物理方向；distance是分配区域较短边一半的比例，以中心面板的分配区域中心为起点。

{{code:courseblockbn}}

![显式毫米布局](../assets/coursephysicalbox.png)

数值box按 xmin、xmax、ymin、ymax给出比例；绝对单位box保持真实物理尺寸，坐标相对于扣除原生外层装饰及连接图例后的内容区。box不与side/size/angle/distance同时使用，center_box不与center_size同时使用。每次目标尺寸均重新检查物理区域，不暗中夹紧。`ggchord_fusion()` 或 `magnify_ggchord()` 的显式 `allow_overlap = TRUE` 允许有意叠放，后加入的面板盖住先加入者。使用矢量grid裁剪，不截取PNG；透明连接不在同一区域重复叠画。导出的 `overlaps` 表记录覆盖顺序。这允许有意遮罩，不代表自动避障。

Sankey来源同时出现在左右两列，需从 `get_chord_layout(plot)$display_instances` 明确选择 `q_instance`/`s_instance`，歧义时不能省略。Radial、Chord使用各自已注册骨架frame。q/s区间、单位和来源记录保持不变。同一显式 `source_hsp_id` 重复显示时，所有参与的Fusion图层（含嵌套内部）均须设置 `allow_repeated = TRUE`；`source_observations` 区分明确来源ID和显示次数，不通过比较表值猜测观测身份。Dotplot有独立双轴，当前尚未注册为Fusion端点。
## 外层 Guide 与面板空间预算 {#outer-guides}

各子图保留独立 Scale 和 Guide。组合层的 Scale 与 `guides()` 描述跨面板连接带；
原生 Guide 可以分别放在不同外侧，无需全部共用同一侧：

{{code:courseblockbo}}

![多侧图例预算](../assets/courseguidefusion.png)

`layout = "fixed"` 保留已有位置；`"adaptive"` 复用原生 Guide 候选/测量引擎，
联合考虑各侧预算和面板区域，再在剩余物理空间重建各子图。可以移动或重排未锁定
Guide；显式位置、方向、行列、字体、key 尺寸和绝对面板框仍是约束。例如宽矮
画布上的默认 Identity 尺可以转为横向，按原来六个标签重新测量长度，不丢掉95或100。

`metadata$guide_layout` 报告选择及锁定项；`guide_bounds_mm`、`guide_sizes_mm`
记录完成图例组合后、测量面板前的最终外层框。嵌套组合的 Guide 报告位于
`diagnostics`。无法满足约束时明确失败。`theme(legend.position = "none")` 沿用
ggplot2 语义隐藏组合图例。组合内部 Guide 暂不支持；各独立子图仍可拥有局部
Guide。这是有界物理空间搜索，不声称全局最优。
## 自动外围排布 {#auto-perimeter}

用 `fusion_panel(..., side = "auto")` 将超过八个独立面板排在中央图外围。
排布在外层 Guide 预算后的实际毫米空间进行，在矩形画布上按输入顺序顺时针
放置。显式尺寸和全部固定方位、角度、box 保持不变。未锁定默认容器可适配，
最低15毫米；这不缩小字体，也不保证任意文字均能容纳。各子图保留自己的标签
及物理尺寸策略。

{{code:courseblockbp}}

![同一 HSP 的十二次有意展示](../assets/courseautofusion.png)

真实示例刻意将一条已观测 HSP 展示十二次，显式 `source_hsp_id` 和
`allow_repeated = TRUE` 区分十二条显示连接与一条来源观测。每个视图仍按实际
尺寸重建精确区间端面。`metadata$panel_layout` 记录有界搜索、保留约束及容器
适配，面板行提供 `auto_angle`、`auto_size_factor`。多个 magnify 细节可以在
命名方位向量中重复使用 `"auto"`，同时明确 `detail_panel`。无法摆放时明确
报错，不丢弃面板。这是确定性保守搜索，不声称全局最优或任意数量都能放进任意画布。
## 原生外层图形装饰 {#plot-furniture}

使用标准 `labs()` 和 plot Theme 元素配置组合标题、副标题、图注、tag、背景及
边距。这些设置此前只进入样式 carrier、未被绘制；现在按实测空间预留并出现在
最终静态图中。各子图仍保留自己的轴、主题和源坐标，Fusion 本身没有笛卡尔数据轴。

{{code:courseblockbq}}

![外层标题与图注](../assets/coursestyledguidefusion.png)

`metadata$content_bounds_mm`、`content_canvas_mm` 记录实际内部组合区域，
`plot_decoration` 保留生效标签及原生物理包围框。box 相对于扣除外层装饰和 Guide
后的内容区；已有固定框在目标尺寸重新检查，不够时扩大画布或调整显式边距。
使用原生 `labs(title = NULL, subtitle = NULL, caption = NULL, tag = NULL)`
清除标题。没有装饰时保留旧渲染路径。

嵌套组合保留独立标题及毫米位置注册。原生前景/背景分离和继承裁剪可避免显式
重叠模式下半透明连接重复叠画；源表、已选区间及颜色映射保持不变。

普通子图边距不会充当组合连接的遮罩。显式重叠模式下，连接在前景下方穿过完整容器，后层容器仍按声明顺序覆盖。需要为标注留空时使用实测 `route = "avoid"`，不依赖意外的边距断带。
