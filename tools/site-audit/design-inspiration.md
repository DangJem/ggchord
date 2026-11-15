# 设计来源核对（2026-11-15，v0.28.0）

本记录核对 README 的 Design inspiration。设计文档能证明当时明确列过比较对象，
不自动证明运行过外部软件、采用其算法或完成集成。源码中的同类功能也不能单独证明
历史上的因果借鉴。下表只按现有可追踪记录描述。

| 对象 | 核对依据 | README 应采用的范围 |
| --- | --- | --- |
| ggplot2 | `DESCRIPTION` 的直接依赖；`R/ggchord-build.R`、`R/scales.R`；v0.9.0 的 grammar 整理 | 直接采用分层绘图生态及 aes/scale/guide/geom/coord/theme 职责 |
| BLAST | `R/data-import.R`、`examples/blastn-standard/`；字段 qaccver/saccver、qstart/qend/sstart/send | 实际读取已有比对结果及借鉴紧凑命名；绘图不执行 BLAST |
| Circos | `design-roadmap/versions/v0.10.0.md` 明确写入环形信息层级和 ribbon；环形布局/连接实现 | 明确的视觉设计参考，不能声称复用了其算法 |
| circlize | 仓库规则及 `design-roadmap/versions/v0.17.0.md` 的轨道比较；`R/track-grammar.R` 等轨道实现 | 多轨道扩展能力与 track/ring 契约比较，不能写成复制其接口 |
| clinker | `design-roadmap/versions/v0.10.0.md` 明确记载发表级基因箭头；gene 几何 | 明确的基因箭头/共线性表达参考 |
| gggenes、genoPlotR | 仓库规则及 `design-roadmap/versions/v0.14.0.md` 的线性比较 | 直线基因簇、链方向和共线性语义的比较对象；不声称分别复用过布局算法 |
| gggenomes | `design-roadmap/versions/v0.14.0.md` 的 2026-10-28 调研逐项区分 sync/align；focus、同步、锚定、翻转和 feature stack 的实现 | 已落实工作流概念参照；保留 ggchord 自己的输入和接口 |
| DNA Features Viewer | `design-roadmap/versions/v0.10.0.md`、`design-roadmap/versions/v0.13.0.md` 明确记载特征可读性/冲突处理 | 注释与局部布局的设计参照，不宣称复制算法 |
| ggmagnify | [长期设计](../long-term-audit/DESIGN.md) 明确记录二维区域机制参考及不复制代码；`R/magnify.R` | 概览—细节关联的表达参考；实际来源区间/跨坐标投影由 ggchord 自己实现 |
| ggmsa | [v0.20.0 调研](../long-term-audit/SEMANTIC-RESEARCH-v0.20.0.md) 当时仅列比较方向；本轮实际查看官方 geom_msa、化学配色表及用户给出的矩阵图 | 旧版只能说调研；v0.28.0 有新的具体落实证据：连续单元格、等比例紧凑矩阵及化学性质配色；用户要求深色粗体和更大的字形占比，这是 ggchord 的可读性调整，不声称逐像素复制 ggmsa。配色来源、上游 blob 和 Artistic-2.0 记录在 [NOTICE](../../inst/NOTICE) |
| ggseqalign | 长期设计官方能力表及 v0.20.0 调研；未发现桥接实现或外部工具重测记录 | 仅能力调研对象，不声称已影响每个布局选择或完成集成 |
| ggseqlogo | 同上，记录为未来 logo 可选对接方向；当前没有 logo 图层 | 仅序列 logo 调研对象，原文把它列为既成设计影响过强 |
| Jalview | 长期设计记录 annotation 页面当时访问超时、待复核；v0.20.0 明确沿任务书列为目标 | 仅调研目标；没有证据声称当时完成官方行为复核或采用其算法 |
| ESPript | 长期设计官方能力表及 v0.20.0 调研；原 README 未列出 | 可补充到调研对象段落，不能因此提升为已实现的设计来源 |
| ggview | 仓库规则、静态指定尺寸预览设计；`view_ggchord()` 实现与标准 ggsave 工作流 | 指定真实尺寸预览的工作流参照；保留标准保存接口 |

结论：原列表不是全部凭空编造，但把“已确认的设计参考”和“任务书/能力调研对象”
放在同一语气下，夸大了部分项目的落实程度，并遗漏调研项 ESPript。本轮按证据修订双语源；没有补写并不存在的历史实验。按用户后续要求，当前设计灵感正文不列商业软件产品名称，原始来源与历史审计按各自记录保存。
