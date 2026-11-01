这条完整工作流使用源码仓库已有的四条序列记录，不需要重新比对或预测基因。
使用匹配的开发 checkout，或[下载教学输入](../assets/tutorial-inputs.zip)，解压后让
工作目录下包含 `examples/`。保存图形时同时记录包版本与[构建清单](../build-manifest.json)。

## 读取并检查原始记录 {#read}
{{code:realimport}}

FASTA 头部提供原样 accession version 与长度，GFF3 提供 CDS 区间和真实链方向，
BLAST 文件提供已有比对；读取器不运行 BLAST。继续前检查 `real_validation`：
连接端点递减可能是合法的反向命中。不能把未知注释链改成 `+` 来满足基因快捷层。

## 明确显示筛选规则 {#selection}
{{code:realprepare}}
{{results:real}}

显示连接保留相似度 ≥80%、已报告比对长度 ≥500 bp 的记录。弦图每条序列选择最长
六条已导入 CDS，并只标记其中最长两条。这是显示规则，不是生物学结论，未显示的
基因不表示不存在。原始 `feature_real`、`link_real` 不变，`real_prepared$report`
和 `real_filtered$report` 保留准备与筛选证据。链方向不重编码；这与生成说明中
明确进行显示链重编码的小型随包教学 fixture 不同。

## 构建弦图比较 {#chord}
{{code:realchord}}
{{figure:realchord}}

先读序列名称，再读真实链方向的 feature 箭头及区间连接。ribbon 表达已有比对区间
和映射的相似度，连接交叉本身不构成重排证据。选定注释与更宽的比对区间分开阅读。
科学解释涉及显示子集之外的注释时，应使用完整输入表。

## 在线性视图检查全部注释 {#linear}
{{code:reallinear}}
{{figure:reallinear}}

这张图绘制全部导入注释。ribbon 默认只连接相邻显示序列，保留所有比较需要显式
图层设置。bp 参照条属于共线性坐标。裁剪、源坐标轴与独立窗口见[局部视图](views.html)，
按图集明确尺寸预览保存见[输出](output.html)。

## 随图保留证据 {#provenance}
教学输入包包含未改写的 FASTA/GFF3/BLAST、SHA256 清单和共享示例脚本。
绘图清单记录代码 hash、PNG 设备、宽高和 dpi。accession 名称标识来源记录，
图片不替代原始注释或比对证据。
