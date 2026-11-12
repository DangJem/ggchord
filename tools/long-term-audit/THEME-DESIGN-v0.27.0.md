# 主题系统审计与方案

基线56f76af；2026-11-12结合用户讨论纲要、当前源码与ggplot2 4.0.3核对。
本文是持续开发中的设计证据，不将待实现项写成已完成。

## 已有事实与需要修正的缺口

- ggchord()已隐式添加theme_ggchord，四个少量预设已存在，共享同一实现；无需新造预设。
  默认背景、层级、物理字号保持稳定。Dark应修复跟随主题缺失的非数据装饰颜色，
  不能为了深色背景任意重映射数据颜色。
- theme_ggchord有大量角色图例便捷参数，轴的物理设置存于属性，不全是注册元素。
  已注册的ggchord.axis.*、标签/引线元素可通过标准theme()修改；尚未注册的物理设置
  不能据“theme可填写”认定生效。先补一致入口，现有便捷参数保留兼容，不草率删除。
- GeomGgchordTrack把数值轴颜色#4F5962、线宽.25、字号2.2、字体和角度写死，
  使Dark及用户text设置不一致。数值轴标题还与刻度共用不可区分的绘图分支。
- Track已经有按ID独立的数值、fill、colour scale及guide；同ID图层共享数值域/物理band，
  不应另造主题域或再次复制scale系统。应使用现有scoped scale示例证明其独立性。

## 职责与优先级

| 组件 | 拥有的配置 | 不应承担的工作 |
| --- | --- | --- |
| Data/Stat | 真实记录选择、统计样本、聚合、逐位点/分段证据 | 不因主题更改数据/筛选/HSP计数 |
| aes/Scale | 数据映射、数据单位/值域、palette、breaks和labels | 不以主题色替换用户映射 |
| Guide | 单Guide位置/方向/行列/顺序/形态、局部装饰主题 | 不改scale的数值边界或顺序 |
| Geom | 骨架形态/单双线/带、Ribbon箭头/渐变模式/透明度、Feature形状与固定美学 | 不把数据颜色伪装成主题配置 |
| Position | Feature堆叠/偏移、track band及局部空间调整 | 不改生物学链或来源坐标 |
| Coord | 序列排列、方向/对齐及投影 | 不承担背景/字体或改变关系选择 |
| Theme | 非数据文字、轴/引线、背景/边距、共享图例家具与视觉层级 | 不启动方向/排序算法、不插入新的scale或改变局部Stat |
| Layout策略 | 按实际设备重新测量标签/Guide/面板，在显式约束内搜索 | 不覆盖明确字号/位置，不宣称无限空间或全局最优 |

固定Geom与映射美学遵循ggplot2各自入口，主题仅填非数据装饰默认。
完整预设重置外观；之后的theme()逐字段合并。局部Guide主题覆盖全局同类装饰，
同Track ID有一个数值轴声明，其余同ID层共享它。显式scale始终优先。

## 原型后拟采用的实现

1. 通过ggplot2原生注册机制补齐有实际消费点的ggchord.axis物理/策略元素，沿用已有
   ggchord.axis命名与绝对unit约定。现有theme_ggchord便捷入口转换为同一主题元素，
   后续标准theme()可覆盖；不把两套配置并列长期维护。角色Guide旧便捷入口保留，
   新文档优先展示guides(aesthetic=guide_legend(theme=theme(...)))。
2. 增加最少的ggchord.track.axis.line/text/title及文字/标题offset元素。
   它们全局继承原有轴/文字语义；按Track局部外观优先复用原生guide_axis(theme=...)，
   不引入按每个Track动态注册的主题元素或新的guide函数。不能正确解释的原生选项明确拒绝。
3. 默认/Minimal/Publication/Plasmid保留既有数据外观；Dark仅修复非数据轴颜色。
   字号/字体改变后沿现有目标设备上下文重新计算。局部Guide只改自身装饰。
4. 补实际render验证：同ID叠加、两个ID不同单位/范围/独立fill guide，跨视图主题切换
   不改源表/路由/统计；显式Geom、scale、guide与后续theme()字段覆盖保持。

这些实现结论仅涉及已确认缺陷和兼容补齐，不据此删除现有公开API。进一步缩短
主题签名属于集中API治理时的迁移决策，必须保留可执行替代与删除版本记录。

## ggplot2依据

原生元素类型/继承与扩展注册：
<https://ggplot2.tidyverse.org/reference/register_theme_elements.html>。
标准theme字段合并：<https://ggplot2.tidyverse.org/reference/theme.html>。
单Guide局部theme、position、direction、nrow/ncol/order：
<https://ggplot2.tidyverse.org/reference/guide_legend.html>。
核对的是公开扩展机制，不复制上游实现或依赖未导出的函数。

## 第一工作包实现结果

主题原型已落地，109形参保持不变，没有删除公开入口。axis短参数与原生注册元素
使用同一消费点；同次调用重复拼写报错，后续theme()覆盖生效。Track按原生GuideAxis
局部theme解析，继承也使用ggplot2元素树；只支持能准确解释的theme/angle/check.overlap，
其他选项明确拒绝。数值标题独立于刻度文字；同ID叠加保留一个轴，不同ID仍独立scale。
未锁定的默认文字/标题偏移按最终字形真实尺寸留出间距，显式offset保持用户值。
竖向Collinear/Sankey的物理Track向量按最终屏幕方向换算，避免宽高比例影响物理尺寸。

check-theme-system.R通过真实pUC19的GC比例/A碱基计数、同ID叠加、两种不同单位/值域、
独立colour/fill图例、五个主题预设和四个多序列视图的直接构建。数据/统计保留，
局部Guide覆盖、全局后续Theme、隐藏和不支持输入有定向证据。
第一张代表PNG已查看；大字标题的间距按实际度量修正，不将这些局部保证扩大为
全局数据图层/标注无冲突。本阶段开发记录不代表v0.27.0或全部长期任务已经完成。

后续Fusion审计已修复另一处无效入口：组合labs标题/图注/tag及原生plot Theme
外框实际生效，字体/背景/边距由Theme、文字内容由labs、面板几何由组合负责。
内容空间和Guide联合预算后才建子图；原生轴样式留在子图，组合无虚构笛卡尔轴。
具体定向证据见DEVELOPMENT第七工作包；保持109主题形参，未引入碎参数或删除API。

实际绘制再确认边距不能充当关系过滤：重叠模式中普通子图边距导致断带的反例
已修复，连接跨完整容器并保持独立前景/覆盖顺序。Theme边距负责外观空间，
需要标注留空时由显式路由及真实文字/Guide障碍负责，不隐式裁断科学关系。
