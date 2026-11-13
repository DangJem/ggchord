{{code:courseblockbw}}

{{prepare:fusiondata,fusionpanels,fusionlinear,fusionview}}

## 嵌套源区域保留独立命名空间 {#nested-detail}

细节图可以本身是Fusion。本例选择未改动的pUC19来源上的两个真实坐标窗口，内部主图是较大窗口，内部细节是较小窗口；外层选区通过 `detail/center` 指向内部主图，并非指向图片中的任意矩形。

{{code:courseblockbx}}

![具有叶路径的嵌套细节](../assets/courseouter.png)

`composition_id` 与 `region_row` 将派生选区、front和引线关联到 `focus_region_data_by_composition` 的各层原始选区。顶层 `focus_region_data` 仍是原始外层输入。叶面板用 `detail/center`、`detail/detail` 等路径命名，组合容器不会被当作单一基因组轴。来源在多个叶面板重复时，应指定完整路径与准确实例；`source_panel` 可选择 `center` 下的源叶面板。导出保留实际毫米边界和原始box/angle请求。

各子图保留独立图层、Scale和统计总体，组合器不会暗中复制图层或重算统计。跨原点区间采用后文的显式完整源长度声明。完整负链CDS→基因组focus→真实字符示例见 `examples/m3-nested.R`。字符保留按物理尺寸的tile/error策略，较小嵌套区域不能仅凭增加DPI容纳60个字母。完整多设备发布验收与更广泛的联合布局优化属于独立工作。
## 显式跨原点选区 {#origin-detail}

圆形排版不是环状来源的证据。在 `start > end` 的选区上设置 `wrap = TRUE`。
两个选定叶面板必须保留相同的完整源长度。字符适配器还需显式 `circular = TRUE`，
并在 `source_data` 提供完整字符串或精确 `source_length`；缺少完整字符串时，字符
仍标为未核对。最多覆盖完整一圈，保证源位置到字符列的反查唯一。负链按自身方向
环绕并保留 gap。不会暗中拼接不同的 focus 实例。

完整真实流程见 `examples/origin-detail.R`：从未改动的随包 pUC19 FASTA 提取
82 个实际碱基（2640–2686，再到 1–35）。`focus_region_data` 保留原来的一行；
派生表以 `region_piece` 区分两段，共享同一 `region_row`。只有外侧两个边界有引线，
接缝不是另一条选区或生物学关系。

{{code:courseblockby}}

![显式跨原点选区](../assets/courseorigindetail.png)
