# 内部状态的有限直接验证

来源：`08f753d5b1061fdbaa9ee4eeb91b41f5c8c6782f` 的 R 源码；本批仅另改
`R/zzz.R` 顶部两行注释，执行逻辑未改。环境为 R 4.6.0、ggplot2 4.0.3。

单一最小脚本使用 pkgload 加载当前命名空间，直接核对：

- 在 outer 测量上下文中进入 inner，inner 报错由调用方捕获后恢复 outer；outer 返回后
  恢复原有 context。核对返回结果与错误保留，没有更改恢复策略。
- 一个只含基础序列层的 plot 经 `ggchord_render_plot_copy()` 复制后，修改副本 ref 的
  哨兵字段、向副本 scale 容器加入手动颜色 scale，原 plot 的 ref 与 scale 容器不变。

结果通过。脚本未绘图或测量任何 panel，没有开启发布设备矩阵，没有运行 testthat、
包检查或 vignette。该结果只证明上述嵌套/报错恢复和 ref/scale 容器隔离，不证明
layer/coord ggproto 深复制、所有异常恢复、所有设备/字体效果或重复 build 的全域纯性。
相应高风险工作由 F26/v0.16.3 的状态与稳定输出阶段继续验证。
