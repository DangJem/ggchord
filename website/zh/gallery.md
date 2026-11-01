## 两种语言共用一个来源 {#source}
六张图均使用[共享文档示例源码](../assets/documentation-examples.R)
中定义的**合成教学输入**。未为改善外观改写任何原始 FASTA/GFF3/BLAST 记录。
两种语言、README 和 vignette 引用相同图片字节，图注分别说明含义。
源码及资产 hash 记录在[资产清单](../assets/plots.csv)中。

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
