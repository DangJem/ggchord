基于 ggplot2 的 R 包，用于**多序列弦图、单序列圆图和直线共线性图**。
使用同一套图层语法，把区间表转成便于阅读的序列比较图。

这里是 **0.16.5 本地开发文档**，基于 0.16.4 候选 API。
当前正式发布版仍是 0.15.0；此处不表示新版本正式发布，网站部署是另行授权的动作。

## 选择视图 {#modes}
{{modes}}

## 绘制第一张比较图 {#quick-start}
先使用显式构造的教学数据。下列标识和数值均为合成示例，不代表生物学观测。

{{code:setup}}
{{code:chord}}

在 R 中打印 `p_chord`，或[按显式尺寸保存](output.html)。依次阅读
[安装](installation.html)、[数据与坐标](data.html)、[三种模式](modes.html)、
[常用调整](adjustments.html)，再进入[进阶专题](tracks.html)和 [FAQ](faq.html)。

## 查找工具 {#discovery}
[函数参考](reference.html)列出全部导出函数、拼写别名及上游重导出。
[迁移指南](migration.html)提供已删除调用的当前替代写法。
[图集](gallery.html)把每张图连接到生成代码和输入。
