# 数据输入和再生成入口

随包对象在 `data/`，内部目录在 `R/sysdata.rda`。读文件示例通常使用仓库的 `examples/` 或安装后的 `system.file("extdata", ..., package = "ggchord")`。重建数据时从仓库根目录运行对应脚本，先核对所需输入。

| 生成器 | 输入与先决条件 | 结果与核对方法 |
| --- | --- | --- |
| [基础教学表](generate_example_data.R) | 仓库内的序列、注释、BLAST 比对及质粒 FASTA | 生成 `seq_data_example`、`gene_data_example`、`ribbon_data_example` 及质粒教学对象。教学注释同时保存展示链和 `source_strand`，后者记录来源方向 |
| [蛋白教学表](generate_protein_data.R) | `examples/protein/` 中完整蛋白 FASTA 和已观测的比对表 | 两个 `protein_*_data_example` 对象，检查每条序列的字符数与全长 |
| [网站数据](generate_website_data.R) | 仓库内的 FASTA、GFF3、比对、质粒与格式样本，详见[文件和 MD5](website-data-sources.tsv) | `genome_comparison_example`、`sequence_detail_example`，完整来源与显示子集可分别访问 |
| [综合案例](generate_flagship_data.R) | 先生成网站数据，再读取内部目录及随包参考对象。[输入摘要](flagship-data-sources.tsv)列出依赖对象 | 每个综合案例的绘图输入独立保存到 `data/`。摘要对应生成时读入的文件 |
| [通用特征目录](generate_common_feature_database.R) | 需要指定版本的本地工作簿、13 份源图谱及跟踪的历史来源 RDS，开发依赖为 readxl、digest、xml2 | 脚本先检查工作簿 SHA-256 和历史记录，再更新内部目录。所需文件名和期望摘要写在生成器中 |
| [引物目录](generate_primer_database.R) | 使用脚本中列明来源的工作簿，开发依赖为 readxl、digest。按脚本固定 SHA-256 核对输入版本 | 更新内部引物目录并保留其他内部对象 |
| [限制位点目录](generate_rebase_database.R) | 本地 `examples/rebase/` 的源目录及等价名称表，普通 clone 需另备输入 | 更新内部目录，校验 4,909 条记录和 4,907 个不同名称 |

`prepare-*` 脚本用于特定来源样本或格式样例，所需外部文件、程序及输出由各脚本声明。使用已有随包表绘图时可直接加载该对象。若要重新生成外部目录，先取得相应来源和版本并满足其使用条件，避免用名称相似的文件替代已核对输入。

早期单序列存档 `single_genome_example`、`single_gene_example` 和 `restriction_site_example` 的原始输入位于 Git 提交 `2939f848` 的 `examples/single-genome/`。
