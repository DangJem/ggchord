这些辅助入口负责检查、选窗和结果消费，不替用户推断生物学关系。每一步保留源表，按明确问题选择工具。

## 准备本页 {#input}

本页独立加载 ggchord 和可选 dplyr；使用随包真实序列、HSP 与 pUC19 GC 窗口。

{{code:comparescaledata}}
{{code:trackdata}}


## 目录不是搜索结果 {#catalog}

primer_catalog 返回寡核苷酸目录；aliases=TRUE 返回来源别名。目录方向是描述提示，实际结合方向由 find_primer_bindings 在目标 DNA 上计算。同序列多名称不应被当成多个不同寡核苷酸。按需查看目录后，去[引物课](primers.html)搜索和显示。

{{code:comparecatalog}}
{{console:comparecatalog}}


## 目标尺寸诊断与明确顺序优化 {#diagnostics}

diagnose_ggchord 检查一个指定 PNG 尺寸，返回测量、文字与 guide 提示，不修改原图。export_ggchord_layout 按层导出数据；original_data=TRUE 保留源行。get_chord_layout 返回完整公共布局。optimize_ggchord_layout 返回建议顺序、方向和前后分数，不改输入也不自动把建议应用到图上。优化的目标是连接可读性，不是生物学上更正确的基因组顺序。

{{code:comparediagnostics}}
{{console:comparediagnostics}}


## 轨道与窗口使用相同源坐标 {#track-focus}

这里聚焦 pUC19 的 500–1500 bp，并按 interval 裁剪同一批实际 GC 窗口。focus_ggchord_track 不重新计算窗口均值，不在边界创造观测，也不会用一条线跨过省略的窗口。需要重新计算更细粒度 GC 时，应回到原 DNA 明确新窗口与分母。

{{code:comparetrackfocus}}
{{console:comparetrackfocus}}


## 完整辅助入口学习位置 {#tools}

| 任务 | 入口与输入/输出 | 主讲位置 |
| --- | --- | --- |
| 文件读取 | read_fasta_lengths/read_gff3/read_blast/read_bed/read_paf/read_genbank；原始文件到标准表与来源报告 | [格式](formats.html)、[比对](alignments.html) |
| 受支持对象转换 | as_ggchord_seq/as_ggchord_feature/as_ggchord_link；明确 ID 与完整长度到标准表 | [适配](formats.html#adapters) |
| 校验与显式修复 | validate_ggchord_data/clean_ggchord_data/prepare_ggchord_data；问题与修复报告不同 | [数据工作流](data.html#validation) |
| 基因组/contig 键 | validate_ggchord_hierarchy；只消费显式成员关系 | [层级](data.html#hierarchy) |
| 连接整理 | filter/deduplicate/merge/bundle_ggchord_ribbons；保留成员、估计标记和报告 | [连接处理](ribbons.html) |
| 目标与窗口 | find_ggchord_loci/focus_ggchord_data/segment_ggchord_data/map_ggchord_position；显式选择与可逆映射 | [真实目标](feature-geometry.html#locus-selection)、[局部窗口](views.html) |
| 同源与编码投影 | link_ggchord_clusters/build_ggchord_cds/project_ggchord_features/project_ggchord_links；需成员、CDS 与相位证据 | [同源](homology.html)、[投影](projection.html) |
| 生物注释搜索 | primer_catalog/find_common_features/find_primer_bindings/find_restriction_sites/filter_restriction_sites；目录、匹配和显示集合各有含义 | [识别](recognition.html)、[引物](primers.html)、[酶切](restriction-sites.html) |
| 预览与输出 | view_ggchord/diagnose_ggchord/export_ggchord_layout/get_chord_layout/optimize_ggchord_layout；按目标尺寸验证与消费结果 | [输出](output.html) |




## 练习与替换输入 {#adapt}

只改变诊断画布宽度，查看测量提示如何变化；再用另一个真实源区间聚焦同一轨道，检查 source 元数据。用自己的数据时保留未修改的原始表，并在每个处理步骤后查看 report。不要按固定列位置消费布局或把绘图中的 layer_id 当作永久生物标识。
