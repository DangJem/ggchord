# v0.26.0 — M2 字符细节与真实 CDS 翻译

状态：约定范围的本地实现、定向验证与相关文档完成，不是正式发布验收。
2026-11-12 按用户明确授权从整合后的本地 main（03805fa）创建 v0.26.0；使用原
OneDrive checkout，没有 worktree 或仓库副本。DESCRIPTION 为 0.26.0.9000，
最近正式版仍 v0.16.6。没有合并 main、push、部署或 GitHub Release。

## 原型、架构和实际交付

先复用普通 ggplot build/gtable 与 ggchord_panel_context，证明单面板 Cartesian 图
可以在 Fusion 容器内按实际尺寸测量；同时以 EMALE BVI_023A_0009 负链真实 CDS
证明现有 spliced map 提取的 DNA 交给 Biostrings 后与上游蛋白逐字符一致。
原型通过后才确定两个公开入口，没有新增 Coord、另一套源坐标系统或翻译引擎。

- `register_ggchord_alignment(plot, data)` 向普通可编辑 ggplot 附加已显示的字符表。
  x 为原 alignment_column、离散 y 为对齐行 accver、label 为真实 symbol。保留 source
  ID/position、gap、strand、unit 和任意来源元数据。Fusion/magnify 在最终实际面板内
  重建并核对注册字符是否仍位于原位置，identity/reverse x 标度可用，显式移动字符、
  非线性/翻转坐标、分面及缺失源区间会拒绝。共享来源的多个行必须指定 detail_instance。
  宏观源区间对应端点真实字符的半格边界，内部 gap 不被插值成碱基。字符 front 只保存
  精确的两个边界，不把采样值伪称逐位点对应。可与既有 Collinear 多独立细节混用。
- `translate_ggchord_cds(cds, sequences, alphabet, initiation, ambiguous)` 复用既有
  build_ggchord_cds 的分段/链向/显式顺序、phase、partial 和 translation-exception
  契约。完整正向 DNA/RNA 字符串按生物学分段顺序提取，负链互补；仅初始 phase
  和不完整末端被排除，后续 phase 校验连续性，跨 exon 的 codon 不丢碱基。
  Biostrings（既有可选依赖，Artistic-2.0）提供翻译表和翻译计算，没有复制其代码。
  每个模型必须有已知翻译表；literal 为默认，start 使用替代起始规则且要求完整5′端；
  ambiguous 默认 error，显式 solve 不能唯一解析时为X。内部/末端终止符均保留为*。
  models 报告原 partial/phase/remainder、编码及蛋白字符串和策略；nucleotides 保留
  基因组位置、codon_index/codon_base、源 CDS 行/段；residues 的列表列逐 codon 保存
  三个真实源位置/行/段。缺失/过短来源、未知表、异常翻译声明均明确失败。
- 修复 geom_residue 首末字符未参与完整格边界训练的问题；现在用半格范围训练位置
  标度。字体、颜色语义与 tile/error 策略不变，外观变化原因已写 NEWS。
- 字符图仍独立拥有映射、图层、scale 与统计。引线只是装饰，不加入关系/HSP表。
  biological Fusion links 也可使用注册端点，但须明确 fill；字符没有骨架颜色供渐变。
  实际位置测量与 draw-time 字形可读性分开报告，不把布局导出当作字形已成功绘出的证据。

## 真实来源与许可

原 DNA/AA gapped 输入沿用 examples/alignment 的 Biopython 官方输入、许可和hash。
默认 source_position 仍是去gap后的输入片段位置，宏观图使用同一片段长度；没有宣称
这是完整基因组。反向互补对照是对同一真实对齐显式重编码，不冒充独立链向发现。

新 examples/cds/NC_000932.gb 是 Biopython 测试档案保存的原始 GenBank NC_000932.1
叶绿体记录（2009历史版本，不声称最新参考），未改写字符或注释。来源URL、取得日期、
Biopython许可和NCBI Molecular Data Usage说明在该目录 README 中。
`data-raw/prepare-cds-inputs.py` 用成熟 SeqIO 提取完整正向 FASTA 和原注释中的
rps16（负链分段）、psbK（正链完整）、petD（正链分段）坐标、明确codon_start、
表11、蛋白翻译及qualifier。源文件hash固定检查；生成物hash和预期蛋白在 expected.json。
没有运行预测、比对或从端点制造序列。临时安装的 Biopython 仅用于数据重建，不修改
用户 Python 环境或 R 运行时依赖。EMALE 的原始输入与来源记录保持不变。

## 定向验证证据

`Rscript tools/long-term-audit/check-m2.R` 通过；结果在 v0.26.0-m2-results.json。

- 三个叶绿体 CDS 共303 codon，两个真实跨 exon codon；正/负链源字符逐位点核对，
  每个 codon 的三个源位置与既有公开投影交叉核对。加上 EMALE 负链 CDS，共四条
  翻译与原记录预期蛋白逐字符相等。比较时仅删除一个源 FASTA 未保存的末端*。
- 真实 DNA 两行原321–380列（120单元、2 gap），AA两行原221–260列（80单元、4 gap），
  均保留真实差异；源位置/列精确往返。AA使用reverse x标度；按独立训练标度/Coord和
  导出面板边界核对字符端点毫米位置。两侧源表保持不变，不增加 biological connections。
- 真实 rps16 的60 bp/20 aa字符+翻译中心与完整叶绿体 Circular主图关联。两个exon源区间
  分别连到同一字符面板；包含内含子的外包区间不能冒充可见编码序列，被明确拒绝。
- 显式反向互补原对齐的source核对与负向端点、共享来源歧义/明确实例、修改注册位置、
  非线性标度、完整/partial、phase不连续、RNA、两种翻译表、替代起始、模糊codon、
  缺失/过短来源、不支持翻译表和已声明翻译例外均有直接检查。
- 一次混合字符/Collinear多细节、显式fill的 biological Fusion link、原单Collinear路径通过。
  标准 ggsave 实际生成 DNA、AA、CDS 三张PNG，均人工查看。尺寸16×5英寸、130dpi，
  字体2.7mm，为这些输入的定向证据，不是全部画布/字体/设备矩阵。
- 窄画布明确error及保留tile路径通过。预期draw-time错误关闭临时设备时，R报告
  `Killing locked device` warning；成功的三张PNG没有依靠此错误探针作为绘出证据。

## 文档与剩余限制

NEWS、两个新增入口及相关 roxygen/Rd、数据契约、双语formats/fusion课程、兼容矩阵、
安装版本、API学习路径与参数图谱同步；examples/m2-detail.R独立执行通过。
formats与fusion课程各在独立R环境执行通过。formats旧Seqinfo合并警告和fusion旧反向
HSP验证警告保留，没有篡改真实输入来消警告。1059参考译段、生成参数/课程/轻量Rmd
源一致性、双语原始R块与锚点通过。暂存差异空白检查除原始GenBank外通过；原文件
ORIGIN尾空格及末尾空行原样保留，hash核对不变，不为通过格式检查篡改源输入。Rmd只同步源码，没有编译。

M2的本版核心字符注册、CDS真实翻译/phase与宏观关联已经实现，未将它们延后。
AA残基使用CDS ID的已观测aa命名空间；bp宏观图若要关联蛋白，应显式CDS投影或使用
已展示的编码碱基/翻译中心图。分段/跨原点选区须显式拆分，不自动猜缺失序列或相位。
未支持完整MSA分析、CIGAR全格式、翻译例外/移码分析、任意facet/Coord、任意inset/
嵌套或全局引线避障。M3收敛、Fusion F2/F3、更多适配器与1.0治理仍是独立长期待办。

仅增加/复用最小公开成功案例留待正式发布授权后执行；本轮未运行testthat、
devtools::test、R CMD build/check、vignette编译或完整PNG/PDF/SVG发布矩阵。
没有重建冻结docs、旧lesson engine或178旧图，旧远端CI未通过新push验证。
