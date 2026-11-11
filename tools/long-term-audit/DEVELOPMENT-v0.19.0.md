# v0.19.0 Composition 本地开发记录

状态：运行时首批完成，文档与最终配图收尾中；尚未整合本版或发布。
用户明确批准将v0.18.0 b038b68快进整合本地main，并从main建立v0.19.0，持续开发至本版完成。
已执行整合/建分支；使用原OneDrive checkout。本版不得据此再整合main、push、部署或发布。
最近正式版仍为v0.16.6。包开发版本0.19.0.9000。

## 运行时

- Fusion为独立组合层，借用真实ggplot载体兼容ggplot2 4的S7加法分派，不引入新Coord。
- Circular/Chord中央图，最多四个横/纵Collinear外围图，上下左右各一个，可含多条序列。
- 每个面板独立复制来源图与图层，隔离缓存；最终容器内实测gtable与panel，先应用用户原生位置标度与最终Coord，再映射全局毫米。
- 三层绘制为各子图背景、共享Ribbon引擎生成的跨面板连接、子图前景；保留显式背景，
  连接不会被整张不透明gtable吞掉，注释仍在连接前景。无栅格截图拼接。
- 标准q/s字段与可选panel/instance选择，歧义拒绝；来源窗口按原accver和source_origin解析。
  原始四端点、数据属性、元数据及source_hsp_id不改写；重复声明身份须显式授权。
- 新组合前缘采用一基闭区间[min-1,max]并保留方向；连续参考曲线复用Track框架，单位置非零宽。
  导出前缘完整曲线、逐点切/法向、控制法向和全局路径；尖头仍保留真实来源前缘。
- 复用已有Bezier与GeomRibbonFlow，默认Identity填充，渐变/subject尖头显式开启。
  组合层使用普通fill/colour/alpha/linewidth/linetype映射与标准标度；面板角色scale独立。
- 新的available_canvas_mm只在组合容器生效，单图默认上下文不变。Guide实际物理测量后保留空间；
  用户显式主题优先，不能容纳则报告。fixed默认，adaptive仅调整中央物理区域并重测局部布局。
- print/grid.draw、标准ggsave及共享view_ggchord通过同一求解器绘制，不缓存旧尺寸图片。
  导出使用临时测量设备，恢复原设备；实际保存按目标设备重新度量。

## 定向证据

`Rscript tools/long-term-audit/check-fusion.R`通过：
真实EMALE PAF全部11条HSP与原始四端点；Circular/Chord、四方向横纵Collinear、
10×8与12×6英寸的逐前缘投影，来源图缓存隔离、单位置、实例A1/B/A2、
显式scale/theme/guide、源空间Track统计与标签可见性、空关系、未知方向省略报告、
重复source_hsp_id拒绝、错误panel/歧义/不支持视图/小画布错误、改名focus窗口来源坐标。
额外保存临时F0/F1/渐变尖头PNG并人工检查；不是发布设备矩阵。
新增一条最小公开行为测试，未执行testthat。未运行包级检查或vignette编译。

## 支持边界

本版交付F0/F1及联合布局的首批基础，不宣称F2/F3或完整字体缩放系统已完成。
不支持外围Sankey/Radial/Dotplot、任意角/八方位、外部ggplot坐标猜测、跨单位直接连接、
部分可见区间的隐式裁切、跨原点环绕解释、全局聚束/连接避障/零碰撞保证。
原生局部图层继续由各自Coord契约决定；Fusion不为不支持的图层组合自动降级。
跨面板关系允许穿过面板内部、交叉或扭转；现有前景标签保持可见，但不声称全局避障。
Guide保留局部独立语义，不自动合并不同scale或静默隐藏；需要去重时用户显式设置show.legend。
完整固定/缩放/自适应字号三模式与高密度自动路由留到后续独立门槛。

追加定向证据：真实蛋白BLAST保留aa单位与源行；原生scale_x_reverse的圆图反射及来源前缘一致，已确认并修复对称范围下遗漏位置标度的锚点错误。最小Fusion公开测试同时覆盖此回归，尚未执行testthat。

## 最终运行时收敛

- 子图采用独立layout/scale与图层、Geom/Stat/Position副本，不序列化整个quosure/数据函数环境；
  大型教学环境不会因每个面板复制而递归复制无关对象。
- 默认Identity色条按最终摆放方位及真实字体测量长度；底部横向Guide实测约117.6×7.9 mm，
  不沿用旧竖向85 mm高度。明确指定的Guide长度保留并在不足时提示，不强行覆盖。
- show.legend采用标准NA/TRUE/FALSE及命名逻辑向量，原生ggplot2处理图例参与规则。
- 组合构造/加法/打印更新last_plot；内部样式与坐标辅助图不会替换用户的最后图形。
- 最终check-fusion.R全部通过，结果与源码指纹见v0.19.0-results.json。
- 已完整重现178张PNG，原有176张与基线字节相同。随后明确隔离98个未变的非Fusion R文件、
  266段未变的原共享示例与Fusion专用实现/准备段；后续只重绘两张Fusion图。
  这一可复核边界见v0.19.0-figure-isolation.json，共享布局/测量/Guide源码仍使全部图像失效。
