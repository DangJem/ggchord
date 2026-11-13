## 准备 {#input}

{{prepare:foundationdata}}

## 命名 A–B–A 路由 {#instances}

三个显示行来自两个生物学序列。`seq_order` 的值标识来源，名称标识实例。
先筛选每层的原始 HSP，再选择路由；`between` 接受两个名称或两个声明序号。
歧义路由会报错。以下两个子集互斥，合起来保留所有 HSP。

{{code:foundationinstances}}

{{figure:foundationvertical}}

同一个 `foundation_plot` 可以使用 `coord_chord()` 或 `coord_collinear()`。
纵向 Collinear 改变全局视图，保留源链和源坐标，局部正法向朝左。显式文字样式优先。
`map_ggchord_position()` 接受实例名称，纵向反查传入最终 y。

`export_ggchord_layout()` 提供 `display_instances`、`instance_routing` 和
`instance_omissions`。显式 `source_hsp_id` 可检测跨层重复显示；只有确实需要重复时才
设置 `allow_repeated=TRUE`。没有 ID 时身份只在层内有效，不能以行值相同认定同一来源。
