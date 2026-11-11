# v0.20.0 适配器输入与语义审计

本版已定原则是“适配器优先、按实际输入逐项评估”，不是一次提供全部格式。
审计仓库examples文件、data-raw重建规则与当前read_*/as_ggchord_*契约，原始数据不修改。

| 输入 | 现有来源与可落地范围 | 本版结论 / 后续门槛 |
| --- | --- | --- |
| BEDPE | 有原始EMALE PAF11行可作明确格式转写；不是真实paired-read文件 | 正式read_bedpe已实现；坐标往返、任意score、raw/extra/行号、未知锚点审计策略已验证。不制造支持度/identity |
| FASTA DNA / AA | examples/fasta四基因组、examples/protein四完整蛋白；有来源和重建记录 | 显式keep_sequence保留真实字符，默认输出不变；focus同步裁剪。字母表及逐位点展示只在隔离原型显式声明 |
| VCF/BCF/SV | 本仓库没有可验真的VCF/BCF样本与预期解析结果；BEDPE不能冒充VCF | 暂不新增VCF解析器。先选成熟可选依赖与真实SNP/Indel/BND/SV样本，明确REF/ALT、END/SVLEN、置信区间、多等位/样本字段及缺失支持度，不从断点猜区间/同源性 |
| bedGraph / BigWig | 有真实序列可计算明确的GC窗口；没有真实coverage/signal原文件 | 现有Track可接标准区间/数值表；GC只代表实际GC，不能冒充测序覆盖度。本版不新增无真实文件验证的专用reader。后续保留零基半开边界、缺失段、单位及local/global统计口径；BigWig优先成熟可选依赖 |
| InterPro/Pfam domain / motif | 有完整蛋白及BLAST，但没有可追溯的InterPro/Pfam分析结果 | 通用Feature可表达用户给定的AA区间；不把BLAST HSP当域或运行域分析制造结果。先取得真实域表，保留模型ID/AA坐标/evalue/来源后评审轻量适配 |
| CIGAR/真实对齐 | EMALE PAF有cg:Z标签，query BVI_023A FASTA存在；subject Cflag_017A完整字符串当前未随包提供。蛋白BLAST为标准12列，无gapped行/BTOP | 不声称“没有CIGAR”，但当前证据不足以核对两侧真实字符/所有错配。后续可按原归档成员哈希加入真实subject与已对齐蛋白片段，再验证双向column/source映射；六端点仍不够 |
| MUMmer/chain/MAF/MCScanX | 未随包提供真实输入 | 保留统一pair/feature契约，不先铺无样本解析器 |
| BAM/CRAM/GFA/ODGI | 未随包提供、且非本版高优先级 | 按v10使用成熟生态按需桥接，不新增分析引擎或完整组装图Coord |

本轮没有下载新业务数据或引入运行时依赖。新增示例均来自现有有来源记录；合成数据仅用于
明确标注的格式/字母表边界，不作为科研观测。格式解析之外的候选均保持未实现状态。
