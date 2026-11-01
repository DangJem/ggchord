<p class="lesson-intro">按要回答的问题选择案例；缩略图用于找到方向，点击可查看完整尺寸。每张图都链接到数据准备、代码和读图说明。</p>

## 按任务找图 {#tasks}
| 总体关系 | 单序列注释 | 局部位置 | 定量与排版 |
| --- | --- | --- | --- |
| [弦图](#realchord) | [质粒](#plasmid) | [CDS 邻域](#homologs) | [GC 轨道](#multitrack) |
| [点连接](#pointlinks) | [注释堆叠](#featurestack) | [聚焦](#focus) | [标签布局](#labelauto) |

## 两种语言共用一个来源 {#source}
所有图均使用随包数据或有来源的真实记录，见[共享文档示例源码](../assets/documentation-examples.R)
。未为改善外观改写任何原始 FASTA/GFF3/BLAST 记录。
两种语言、README 和 vignette 引用相同图片字节，图注分别说明含义。
源码及资产 hash 记录在[资产清单](../assets/plots.csv)中。

各专题说明显示筛选。GC 来自实际 DNA 计算，同源分组例子核对实际编码序列；
只有诊断示例会故意修改原始记录的副本，并明确标注，不把错误副本当作生物学数据。

## 图片与尺寸 {#figures}
{{gallery}}

## 重建 {#reproduce}
在源码仓库中准备好建议的开发依赖后运行：

```sh
Rscript tools/documentation-plots.R
python3 tools/sync-tutorials.py --check
python3 tools/build-website.py
python3 tools/check-documentation.py
```

PNG 设备仅用于更新这些文档图片，不是输出设备矩阵或发布验收。
源码和尺寸位于 `tools/documentation-plots.R`，参数变化后须重建清单。
自己的图片导出方式见[输出](output.html)。

## 资产来源 {#assets}
[原画、模板和依赖许可](../assets/ASSET-SOURCES.md)。
