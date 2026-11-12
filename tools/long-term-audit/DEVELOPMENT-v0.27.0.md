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
