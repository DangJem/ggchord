# Track 当前语法与统计边界

本次回查v10 §3.6–3.7和实际track-grammar.R，不按已知类名数量宣称全生态兼容。
早期14几何及identity/count/bin/summary有既有定向证据；本阶段在其基础上补齐
smooth/density/ecdf源空间原型、实现及估计组来源，公开签名不增加新碎参数。

| Stat | 可用Geom | Data Position | 源空间与证据 |
| --- | --- | --- | --- |
| identity/count/bin/summary | 原14几何中原生ggplot2输入契约合法的组合 | 已有identity/stack/dodge等原生兼容；bin须显式边界 | 既有validate-track-grammar.R及v0.23严重缺陷回归；不是98组合全认证 |
| smooth | line/path/point/ribbon | identity字符串或原生PositionIdentity | 原始x/y与明确method/formula，ribbon需SE范围；原生lm结果/源组及CI实图核对 |
| density | line/path/point/area | identity字符串或原生PositionIdentity | 事件源位置x，可选weight，带宽由用户声明；真实G/C位置与原生density对照 |
| ecdf | line/path/point/step | identity字符串或原生PositionIdentity | 事件源位置x，累计比例；pad=FALSE显式有限锚点；真实G/C与原生结果对照 |

新Stat其他组合明确拒绝；area默认stack，density/area须主动选择identity。
局部x/y不是最终画布x/y；Stat先于物理band/法线投影，也先于Track数值标度变换。
若需要log空间模型，须准备并命名该变换后的y，不能把后置显示scale伪称先变换拟合。
bounds仍属于源几何error/clip，不兼任原生density的边界校正参数；不静默双义。

来源wrapper通过ggproto方法保留原生参数发现，捕获native compute_group收到的
有效输入行。source_rows是提交给估计器的来源组，不是每点权重/影响、每个拟合点
的新观测数，或自定义模型进一步model.frame过滤的认证。count/bin/summary继续按
各自统计口径保留来源，不能把模型采样点和实验观测求和。缺失源数据原表不删。

## 已有14Geom、支持层次与未公开能力

14类为point/line/path/step/col/bar/tile/rect/area/ribbon/segment/linerange/errorbar/text。
源输入字段遵循已实现局部x/y及范围映射；同ID共享band/值域、颜色标度按ID隔离，
图层绘制次序不等于band order，数据position与physical track相互独立。
本阶段不改变这些已发布/开发契约，也不删除四个高频快捷函数。

- boxplot/violin：原生语义包括类别槽宽、分位数/密度横向宽度及离群点附属几何；
  当前通用路径渲染不是其完整投影契约，不能只把类名加到枚举。已有summary+
  ribbon/errorbar可明确展示用户给定的分位数范围，但不冒充原生箱线/小提琴。
  它们是任务书要求评估的非笛卡尔条件项，本轮不新增不完整封装。
- facet/全局coord_flip/raster：当前来源注册和物理band仅单panel，每个facet都需
  独立source实例/尺度/设备区域，不能透传原生facet后共享错误的唯一参考架。
  独立ggchord子图/Fusion是现有明确可用路径，不宣称facets已经实现。
- Track group为可选高级层级。现有ID叠加、series group、order和scoped scale已经
  表达真实GC/A单位轨道及CI/mean叠加；尚无2–3个需要额外group容器的真实场景，
  维持候选，不增加强制用户层级或为每ID复制一套Theme参数。
- 任意第三方Geom/Stat对象没有投影/来源证明，仍不作为全透传承诺。

## 实际检查

check-track-estimates.R逐Stat与同一ggplot2 4.0.3实时计算对照x/y/CI/density/ecdf，
用真实pUC19 GC窗口与G/C位点分组；不是锁定跨版本固定浮点常量。缺失值格式
控制、全估计组来源、显式有限ECDF、数值scale后变换及Circular/竖向Collinear
保持源统计一致。CI带+均值线PNG实际查看，line/path/step原误取ymax的问题修复。
结果JSON及双语tracks课程独立执行通过。复用最小Track测试只补一条公开smooth
成功路径，没有运行testthat或包级发布流程，未做全输出设备/平台认证。

后续跨视图/布局/全站证据继续，以上不替代v10全部任务完成或稳定API冻结。
