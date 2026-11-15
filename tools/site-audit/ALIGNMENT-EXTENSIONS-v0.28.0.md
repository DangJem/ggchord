# 字符视图附加能力评估（待确认）

2026-11-15 用户询问 ggmsa 的 seqlogo、GC、seed、msaBar 和 helix 是否值得借鉴。
本记录是建议，不是 v0.28.0 的实现承诺；本轮已授权新增的是残基换行。

| 能力 | 官方资料核对 | ggchord 当前可复用内容与建议 |
| --- | --- | --- |
| 序列 logo | [geom_seqlogo](https://yulab-smu.top/ggmsa/reference/geom_seqlogo.html) 展示逐列堆叠字母；其成熟字体和 logo 生态值得对接 | 有价值，适合放在字符矩阵上方，但优先桥接 ggseqlogo 等成熟组件。先确认频率/概率/信息量单位、gap/含糊字符、序列权重和小样本处理；统计结果与几何分开，不从 HSP 或一个 CDS 生成“群体保守性” |
| GC | [geom_GC 官方源码](https://github.com/YuLab-SMU/ggmsa/blob/master/R/geom_GC.R) 在对齐行旁以点大小表示该行 GC 比例，并非滑窗曲线；源码的分母还包含其他输入字符 | 复用已有 `geom_track()`、`geom_track_col()` 和 GC 教程；行汇总可用普通 ggplot 点/条。优先补统一数据整理和字符旁组合范例，避免重复新建同义 Geom。须明确全长/当前窗口和 gap/IUPAC 分母 |
| seed/motif | [geom_seed](https://yulab-smu.top/ggmsa/reference/geom_seed.html) 用底色或星号突出 miRNA seed；[源码](https://github.com/YuLab-SMU/ggmsa/blob/master/R/geom_seed.R) 定位首行首次匹配并高亮对应列 | 高优先级建议：通用 motif/选区标记可用于引物、结合位点和 CDS 局部细节。先复用 feature/来源映射，字符行补明确的对齐列或源区间注释。分开“查找匹配”和“绘制高亮”，不得默认把首行命中位置套到所有序列 |
| 列统计条 | [geom_msaBar](https://yulab-smu.top/ggmsa/reference/geom_msaBar.html) 称保守性条；[源码](https://github.com/YuLab-SMU/ggmsa/blob/master/R/geom_msaBar.R) 实际使用每列最多字符的出现次数 | 高优先级建议：先提供可追踪的列频率/一致率数据，再用普通 `geom_col()` 或现有 track 绘制。次数、有效分母、一致率和熵必须明确，不能统称 conservation 而不定义；序列筛选与换行不得偷偷改变统计口径 |
| RNA 二级结构弧 | [geom_helix](https://yulab-smu.top/ggmsa/reference/geom_helix.html) 接收外部二级结构，能比较已知/预测结构 | 有价值但后置，需要成对碱基/结构来源、配对到对齐列、gap 和伪结处理，以及与现有 HSP 连接的独立语义。接受外部结构结果，不从普通序列或 HSP 自动预测结构 |

建议顺序：通用区间/motif 标记与列统计 → logo 桥接及 GC 行汇总范例 → RNA 配对结构。
是否进入某一后续版本、公开名称和统计策略仍待用户确认。不按第三方函数表整套复制，
不将第三方品牌写入 ggchord 函数或参数名；本轮仅阅读资料，没有复用这些函数代码。
