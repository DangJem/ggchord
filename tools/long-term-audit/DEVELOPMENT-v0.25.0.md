# v0.25.0 — 多独立源区间细节面板

状态：约定范围本地开发完成，未整合、推送或正式发布。
按用户批准将v0.24.0 074899e合入本地main，从该main建立v0.25.0。
main保持074899e，origin/main保持69556d8；本轮没有远端操作。
版本0.25.0.9000新增现有公开接口的多面板能力，没有再拆分额外版本。

## 任务书与实现

依据v10 §4.10.6 M3，在已完成M1与可选尺寸策略的基础上，实现非嵌套、四侧显式多细节。
`magnify_ggchord()`的detail_plot接受按唯一面板ID命名的一至四个Collinear图的列表；
region的detail_panel逐行指定目标，side为按面板命名的不同侧边向量，detail_size可为
命名尺寸列表。单个图的原有参数、默认分配和调用保持兼容。

每个细节图独立持有数据、映射、统计和标签。沿用Fusion区域分配、重叠/超界诊断、
实际毫米测量及统一源区间投影。每个面板至少有一条选区，不猜测缺失目标。
focus_regions增加source_panel/detail_panel；region_row继续连接派生表和原始行。
focus_region_data保留原表及任意元数据。引线仍为装饰，不进入生物学connections或
relationship_data。未新增Coord、公开函数、依赖或数据模型。

## 定向验证

命令：`Rscript tools/long-term-audit/check-magnify-multiple.R`。
使用随包四序列数据及第一条真实HSP的两端，独立聚焦为两个竖向Collinear细节面板，
其中一个施加原生反射标度。检查如下：

- 4×4、8×8、12×8英寸均完成共享构建；逐源行核对两个面板的闭边界。
- 引线目标对应正确面板；独立利用实际骨架、原生标度和实测panel边界核对物理锚点。
  采用相对高亮宽度量级的附着容差，不冻结设备浮点坐标。
- 两个源行产生四条装饰引线，不增加生物学连接；原始同名region_row和前缀元数据保留。
- 同一真实源区间出现在两个独立目标时，按显式panel名称分别选择；不猜测第一个目标。
- 缺少目标、未知目标、重复侧边、显式重叠和超大区域分配均拒绝；旧单图调用仍产生两条引线。

结果见[v0.25.0-multiple-details.json](v0.25.0-multiple-details.json)。
代表12×8英寸PNG已人工查看：[多独立细节](v0.25.0-multiple-details.png)。
引线在该图中会相互交叉，并可能穿过主图；没有宣称全局避障。
复用既有最小Fusion公开案例改为双目标成功路径，未运行testthat。

## 文档

NEWS、roxygen/Rd、英语/汉语Fusion课程、支持矩阵、安装版本及参数图谱同步。
`Rscript tools/check-documentation-lessons.R --page=fusion`在独立R进程通过，
证据[v0.25.0-fusion-lesson.json](v0.25.0-fusion-lesson.json)。
旧课程真实反向HSP的输入验证警告保留，未改动原始坐标来消除警告。
1045参考译段校验通过；62对页面原始R块/锚点一致；参数图谱、课程、轻量Rmd源一致，
git diff --check通过。没有重建网站或旧图。

## 剩余与发布边界

本版只完成M3中显式非嵌套多细节部分。任意inset位置、嵌套、全局引线/标签避障、
自动图层重建/统计口径选择仍未实现。M2字符面板注册和CDS/phase、F2/F3高级Fusion、
1.0前API治理继续保留。完整M3的SVG/PDF文字/mask兼容与发布设备矩阵尚未验收。

没有运行testthat、R CMD build/check、vignette编译或发布设备矩阵。先前远端CI及旧
lesson engine/178旧图指纹事项仍待后续处理，本轮不将其写成已通过。
