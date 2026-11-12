# v0.27.0 持续开发记录

2026-11-12用户明确批准将v0.26.0合入本地main，并持续推进长期任务书所有剩余工作。
已快进main到56f76af，从main创建v0.27.0。原OneDrive checkout，无额外worktree。
最近正式版仍v0.16.6，本轮不推送/部署/发布，不把本地完成当正式发布。

## 主题与多Track第一工作包

已读取用户完整主题讨论纲要，并将12方面约束记录到本版完整设计文档。
源码审计、责任划分、ggplot2公开机制来源和实现取舍见THEME-DESIGN-v0.27.0.md。
默认主题和四个预设已有，复用而非另造；109个形参未新增或删除。

实现序列轴物理配置的标准theme()入口与旧短参数共用、重复入口拒绝；Track数值轴
line/text/title继承与局部原生guide_axis主题，按字形实测留出未锁定默认间距。
支持局部angle/check.overlap；不能解释的原生选项明确失败。全局Theme、局部Guide、
数值Scale、数据颜色Scale与Position band各自独立；修复深色轴装饰和硬编码字体。
同步实际预留空间使用的字体、字号、角度与对齐，保留显式offset。
竖向/镜像物理Track向量用实际显示方向换算。

定向脚本check-theme-system.R使用真实pUC19计算GC比例与A碱基数（不伪造coverage），
验证同ID叠加/不同ID独立scale及单位、局部Guide与后续全局theme覆盖、隐藏、五种预设、
四个多序列视图、源输入与统计保留、错误配置。代表PNG实际保存并查看。
证据为v0.27.0-theme-results.json及v0.27.0-theme-tracks.png。
样式课程独立执行通过；roxygen/Rd、NEWS、两语样式课/安装/参数图谱已同步。
1061译段和两语代码/锚点及生成源一致性通过。预设对数据颜色不作暗中重映射，
深色背景的数据对比度仍由用户geom/scale负责。

这仅是第一工作包，v0.27.0继续开发，不宣称全任务书完成。M3、Fusion F2/F3、
数据适配、全局质量、整站旧资产重现及1.0治理仍在持续范围内。
未运行testthat/devtools::test、R CMD build/check、vignette编译或发布设备矩阵；
冻结docs没有重建，本轮没有远端操作。

## Fusion/M3显式组合第二工作包

先以独立原型核对Sankey/Radial骨架frame与公开map_ggchord_position，再接入共享
Fusion源区间解析。新增八个命名侧边、实际物理angle/distance、相对/绝对单位box、
center_box和显式allow_overlap。角度以中心分配区域为起点，按真实毫米空间定义，
不因矩形宽高比变成另一个角度，也不旋转子图Coord/字体。旧四侧默认保持。

嵌套Fusion递归按实际容器尺寸重建，叶面板用完整路径注册。组合容器不是一个
生物学坐标域；原始源/对齐位置、实例、origin及单位不变。magnify支持detail_box、
center_box、嵌套目标路径和可选source_panel。任意元数据及各层原始选区由
focus_region_data_by_composition保留，composition_id+region_row关联派生表。

显式重叠仍默认禁止。允许后按叶图层顺序覆盖，以互不相交的原生grid矩形剪裁
连接区域，不栅格化子图，不用重复半透明绘制掩盖背景，也不把重叠许可当自动避障。
父子层共享完整渲染账本，保留内部图例、连接和高亮层级。同一明确source_hsp_id跨
嵌套重复显示，全部参与层均需allow_repeated=TRUE；source_observations只按显式ID
区分来源观测与显示次数，不猜同值表身份。

`check-fusion-advanced.R`通过：

- 真实随包HSP在八个独立面板重复展示，Sankey/Radial/横纵Collinear均有代表；
  每个新视图的front独立用公开源位置映射及最终实测panel范围核对（0.02mm容差），
  原q/s端点、来源表及pident保留，八条显示连接仍是一条明确来源观测。
- Sankey两侧歧义拒绝；矩形画布31度物理位置、绝对毫米box、超界/冲突配置检查；
  旧source-region流程在4×4、8×8、12×8定向构建保留源边界。
- 真实EMALE负链CDS→原始基因组focus→60个真实字符，两层原始选区和四条装饰
  引线保留；完整叶路径准确，跨层重复HSP未全部显式允许时拒绝。
- 显式覆盖输出重叠表，默认同区域拒绝；标准ggsave生成八面板和真实嵌套两张PNG，
  均人工查看。布局中的交叉仍可存在，这不代表全局路由/避障已完成。
- 相关Fusion课程在独立R环境执行通过（保留旧反向HSP警告）；新增m3-nested.R
  可运行真实来源示例，帮助/Rd、NEWS、双语课程/参考/兼容矩阵/参数图谱同步。
  最小公开Fusion案例复用为混合嵌套双细节路径，不运行testthat。

证据：v0.27.0-fusion-advanced-results.json、v0.27.0-fusion-eight.png、
v0.27.0-m3-nested.png、v0.27.0-fusion-lesson.json。

第二工作包只完成显式多视图/摆放/嵌套核心，完整F2/F3/M3尚未整体关闭：全局路由/
标签和Guide联合优化、跨原点选择与更多适配器等继续开发，发布矩阵仍未运行。

## 第三工作包：跨原点

新增显式circular字符源转换及wrap选区，默认线性行为不变。来源长度、双向映射、
负链互补/gap和字符核对分开记录，不从图形或HSP推断完整源。原选区保留，
两个绘制片段用region_piece区分，仅外侧起止两条引线。

check-origin-detail.R通过真实pUC19 2640–2686/1–35共82碱基、负链含gap、反查、
缺失/冲突长度与重复绕圈拒绝、精确整数上界及不同实例不可暗中拼接检查。
v0.27.0-origin-detail.png实际保存并查看；v0.27.0-origin-results.json保留结果。
公开示例examples/origin-detail.R可独立运行；帮助、NEWS与两语课程同步。

## 第四工作包：实测障碍与聚束

check-fusion-routing.R通过三条真实HSP的显式共享显示通道、完整front和出发法线、
源记录/观测数保留、完整多边形避障、装饰引线、阻断失败及不同物理画布重建。
原生字体/旋转/图例/标题包围框按目标设备测量，曲线字符串采用保守包络；公开
route/route_width不进入Theme。内部复用ggplot2 compute_just与原生文字保持一致，
这一窄内部依赖需发布兼容验证，不复制第三方代码。

连接与引线导出实际状态、碰撞、通道/成员及折线；嵌套祖先追加覆盖后重新核对，
不能保持clear状态掩盖遮挡。默认direct保持旧几何，不自动过滤。证据为
v0.27.0-routing-results.json及v0.27.0-fusion-routing.png；可运行真实示例
examples/fusion-routing.R。尚不声称全局最优、任意自定义grob逐字形障碍或全发布
设备质量认证；外层Guide联合优化和更广泛密集制图质量仍继续开发。

路由还按实际边框线宽扩张障碍；subject尖头与渲染器共享同一几何函数后核对实际
可见多边形，渐变/尖头和有边框路线已实际绘制。代表图使用竖排真实accession
标签，避免为了展示路由而保留已见的标签相互重叠；原始ID与科学来源不改。

## 第五工作包：外层原生Guide联合预算

fixed保留旧位置，adaptive复用原生Guide候选/训练/测量引擎，加上各侧共同预算、
自动面板尺寸及显式box可行性约束。多原生Guide可分侧，先扣各侧最终框，再在
剩余物理区域重建所有子图。显式位置/方向/行列、字体和绝对box不移动/缩小，
所有源记录及映射保留。默认Identity尺的自动方向和尺寸不当用户锁，转向后重新
按六标签间距计算长度。复用calc_element解析原生主题方向间距继承。修复显式
legend.position=none进入候选失败的问题，沿用原生隐藏语义。

check-fusion-guides.R通过真实Identity宽矮转向（六keys）、双侧原生图例/行列锁、
同一真实HSP表映射不变、面板/Guide不重叠、固定方向/尺寸不够时失败、绝对框不
暗移、inside拒绝、隐藏图例，以及共享单面板自适应引擎代表构建。PNG实际保存，
v0.27.0-fusion-guide-results.json记录结果。metadata最终Guide尺寸/边界独立于
候选carrier测量位置，嵌套diagnostics保留子组合报告。无全局最优或正式设备认证
承诺；自动环绕布局、高密度质量广覆盖、适配器及稳定化继续开发。

## 第六工作包：自动环绕外围面板

side=auto显式选择按实际毫米画布/外层Guide预算的顺时针外围搜索，可超过八图。
固定box/方位/角度和显式size均不动；只有默认容器在有限范围内适配，最低15mm，
不暗改字体、Coord方向、数据或图层。摆放失败明确报错，不丢面板。导出自动角度、
实际容器适配系数及有界搜索策略；嵌套诊断保留子组合panel_layout。

check-fusion-auto.R通过12外围+中央、真实HSP原记录不变/12显示仍1来源观测、
矩形物理方位、框无重叠、混合显式框/尺寸保留、另一宽高比重建、多细节重复auto
方位以及超大显式尺寸失败。真实完整示例examples/fusion-auto.R；
v0.27.0-fusion-auto-results.json与实际保存的PNG为证据。并非全局最优或任意
数量/文字在任意画布的可读性保证；密集质量广覆盖及其他长期工作继续。
