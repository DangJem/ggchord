换坐标系时，先确认数据单位、图层和参数的边界。下表依据当前 Rd 和源码契约；“支持”表示该类图层有对应实现，不表示所有参数组合都经过穷举验证。对照运行记录见版本审计，不把静态说明当成设备验收。

## 图层与三类坐标 {#layers}

| 图层 | 弦图 | 共线性图 | 单序列圆图 | 证据与条件 |
| --- | --- | --- | --- | --- |
| `geom_seq()` | 支持 | 条件支持 | 支持 | [骨架参考](reference/geom_seq.html)；共线性拒绝显式 radius、curvature、ring |
| `geom_seq_label()` | 支持 | 条件支持 | 支持 | [参考](reference/geom_seq_label.html)；共线性不接受显式径向 label radius |
| `geom_gene()` | 条件支持 | 条件支持 | 条件支持 | [参考](reference/geom_gene.html)；要求 bp 与已知正负链，aa 用通用特征 |
| `geom_feature()` | 支持 | 支持 | 支持 | [参考](reference/geom_feature.html)；数据须匹配当前序列与单位 |
| gene/feature 固定标签 | 支持 | 条件支持 | 支持 | [特征标签](reference/geom_feature_label.html)；共线性径向/切向含义改为竖排/横排 |
| gene/feature 自动标签 | 支持 | 条件支持 | 支持 | [自动标签](reference/geom_feature_label_repel.html)；共线性用行布局，拒绝 external=FALSE |
| `geom_link_ribbon()` | 支持 | 条件支持 | 不支持 | [参考](reference/geom_link_ribbon.html)；必须有合法连接；共线性按显示行判断邻接，分段同行带省略并报告 |
| `geom_link_line()` | 支持 | 条件支持 | 不支持 | [参考](reference/geom_link_line.html)；显式分段同行点连接不支持；单圆图拒绝 ribbon/link 图层 |
| `geom_seq_region()` | 支持 | 支持 | 支持 | [参考](reference/geom_seq_region.html)；区间、ID 和侧向须合法 |
| `geom_seq_scale()` | 不支持 | 支持 | 不支持 | [参考](reference/geom_seq_scale.html)；直线标尺只适用共线性图；其他模式用局部序列轴 |
| `geom_track_bar/line/tile()` | 支持 | 支持 | 支持 | [轨道参考](reference/geom_track_tile.html)；必须显式 position_track，同 ID 声明一致 |
| `geom_seq_break()` | 不支持 | 条件支持 | 不支持 | [参考](reference/geom_seq_break.html)；需可识别的窗口边界；不是修改源序列 |
| `geom_seq_center_label()` | 条件支持¹ | 不支持 | 支持 | [参考](reference/geom_seq_center_label.html)；单条序列；正式教学使用圆图 |
| `geom_primer()` / 自动引物标签 | 条件支持¹ | 不支持 | 条件支持 | [参考](reference/geom_primer.html)；圆图与核酸单位 |
| `geom_restriction_site()` | 条件支持¹ | 不支持 | 条件支持 | [参考](reference/geom_restriction_site.html)；圆图、核酸、有效切点 |

本表不将 `geom_feature_plasmid()` 单独当成新几何类型：它是[通用特征的预设](reference/geom_feature_plasmid.html)，主要为圆形质粒教学服务。其他坐标下不要假定其圆图布局参数可原样照搬，先使用通用 feature 图层。

¹ 定向运行发现，单条 pUC19 配合 `coord_chord()` 的中心标签、引物和酶切图层可以 build；现有共线性帮助却将它们称为 circular-only。表中保留这个契约与实现的差异，不把一次 build 当作视觉正确或官方跨模式承诺。引物自动标签尚未单独运行；正式教学仍使用圆图。此任务不修改运行时。

## 读懂条件支持 {#conditions}

`R/coord-collinear.R` 的 `ggchord_validate_coord_contract()` 显式检查圆图参数、自动方向与分段冲突。通用数据契约和核酸专属限制由运行时校验负责。完整源码随网站放在 `assets/source/`；函数参考从当前 roxygen/Rd 生成。

Position 只控制局部放置，不把不支持的图层转换成另一种坐标图层。共线性图的 outside/inside 是上方/下方，圆图是局部法线的外/内；复制配置前应重新读图。显示翻转与基因生物学链的区别见[四种方向](collinear-layout.html#orientation)。

## 版本与复现 {#versions}

本表面向开发版 0.16.7.9000；[安装页](installation.html)说明本地源码安装与正式版入口。单圆图只接受单条显示序列；圆形蛋白只是显示布局，不声明蛋白分子为圆形。混合 bp/aa 数据必须先投影或分别作图。
