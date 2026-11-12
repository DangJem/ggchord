# v0.22.0 — 显式源区间的语义细节组合

状态：约定范围的本地开发与定向验证完成，尚未整合main、推送或发布。
2026-11-12按用户授权将v0.21.0 69556d8快进整合main并成功推送origin/main，
随后从main创建v0.22.0。包开发版本0.22.0.9000，最近正式版仍为v0.16.6。

## 实现与边界

公开magnify_ggchord(plot, region, detail_plot, ...)：Circular/Chord源图与显式
横/纵Collinear细节图组合，四方向位置、独立图层和原生标度，沿用普通Fusion与ggsave。
region是原始accver/start/end正整数闭区间，支持instance/detail_instance消除显示实例歧义。
两个面板须完整包含该源区间，且同为bp或aa。任意图层和统计不自动复制或重算。

共享区间解析、源origin和原生scale→Coord→实测panel毫米投影；高亮宽度按毫米，
两条装饰引线连接源区间的真实边界。高亮在面板前景之后绘制以避免被不透明骨架遮住，
引线在前景之前。普通Fusion没有focus时保持原有绘制顺序。
export_ggchord_fusion新增focus_region_data、focus_regions、focus_fronts、leaders；
原始选择及任意元数据独立保留，派生报告不覆盖原字段。装饰不进入真实HSP关系表。
报告给出源/实例/单位、采样物理长度比和explicit_detail_plot统计口径；长度比可小于1。

本版不包括跨原点单选区、任意画布inset、多级嵌套、多个选区、自动统计口径、
Cartesian字符图接入或全局避障。M1–M3完整目标、Fusion F2/F3、通用Guide/Adaptive仍未完成。

## 定向验证与图像

check-magnify.R通过真实pETDuet-1的101–1900区间（1800 bp），10条feature、2条primer、
3条限制识别区间；多尺寸源坐标保持、反向原生标度、重命名窗口、四源Chord、
重复实例选择、原生aa、单位置闭区间、真实HSP与装饰共存，以及代表性非法输入拒绝。
输入源数据保持；选区原表精确保留。v0.22.0-magnify-results.json、三份CSV和PNG记录结果。
已人工查看最终PNG：选区高亮和两条引线可见；密集细节标签仍可能与局部几何重叠，
不声明全局零碰撞。该PNG是局部功能证据，不是完整发布设备矩阵。

共享解析重构后check-fusion.R定向脚本通过，结果单独保存在v0.22.0-fusion-results.json；
v0.19.0历史结果未改写。已有Fusion最小公开测试补入一次新接口成功路径，未运行testthat。

## 文档

NEWS、roxygen/Rd、汉语参考、两语Fusion教程、参数/任务/API索引、安装版本及Rmd源已同步。
146导出、1604形参、110 Rd；1037段参考翻译校验通过。
Fusion教程在独立R环境执行通过；62对源页面锚点与原始可执行代码一致，124页上下文lint通过。
这些不构成整站构建/浏览器检查；没有编译vignette或重建冻结docs、旧图和网站。

## 推送后的远端CI与本地修复

main 69556d8自动触发的R-CMD-check和Documentation均失败，明细见CI-69556d8.json。
Linux release日志显示197通过、1失败：显式straight限制位点引线被共享排版改成折线。
本地提交a3ba85f保留显式直线，直接重现验证straight为2顶点、elbow为3顶点、内部标签保留；
同时补齐已报告的utils/head、grDevices/as.raster导入与NSE声明。原回归测试保持不变。
修复尚未推送，没有远端重跑，因此不能宣称CI已通过。
Documentation的5个旧lesson engine及178图依赖指纹过期仍待整站重现；部署任务未运行。

未运行本地完整testthat、R CMD build/check、vignette编译或PNG/PDF/SVG发布矩阵。
本版完成仅指上述明确范围；完整v10状态见V10-STATUS.md。
