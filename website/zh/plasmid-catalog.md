这 13 张图使用各自的仓库 FASTA 参考序列。载体用途用于选例，图中注释来自相应的参考序列。完整来源和长度由随包数据及生成规则保留。

## 一次准备全部案例 {#prepare}
共享脚本使用 dplyr 选择序列列。先运行 `library(dplyr)`。

{{recipe:catalog}}

将脚本放入工作目录并运行，然后选一张图跟做。保留全部特征。名称先在圈内放置，放不下的可移到外圈。外圈单独标注最接近三个等距目标位置的实际酶切位点，显示酶名与位置。

## 输入与识别结果 {#counts}
{{plasmidcatalog}}

## pSB1C3：小载体与基本骨架 {#catalogpsb}
```r
sequence <- plasmid_example_pSB1C3
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogpsb}}

## pUC19：常见注释与复制起点 {#catalogpuc}
```r
sequence <- plasmid_example_pUC19
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogpuc}}

## pBR322：多个抗性区与方向 {#catalogpbr}
```r
sequence <- plasmid_example_pBR322
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogpbr}}

## pBluescript II SK(+)：局部调控与短特征 {#catalogblue}
```r
sequence <- plasmid_example_pBluescript_II_SK_plus
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogblue}}

## pET-28a(+)：表达区与标签 {#catalogpet}
```r
sequence <- plasmid_example_pET_28a_plus
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogpet}}

## pETDuet-1：双表达单元，进入密集案例 {#catalogduet}
```r
sequence <- plasmid_example_pETDuet_1
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogduet}}

## pDONR221：重组相关结构 {#catalogdonr}
```r
sequence <- plasmid_example_pDONR221
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogdonr}}

## pcDNA3.1(+)：表达载体的注释层级 {#catalogpcdna}
```r
sequence <- plasmid_example_pcDNA3_1_plus
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogpcdna}}

## pTRE-Tight-BI：双向调控布局 {#catalogtre}
```r
sequence <- plasmid_example_pTRE_Tight_BI
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogtre}}

## pTRIPZ：更长序列的显示密度 {#catalogtripz}
```r
sequence <- plasmid_example_pTRIPZ
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogtripz}}

## pCAMBIA1300：较大植物载体的结构 {#catalogcambia}
```r
sequence <- plasmid_example_pCAMBIA1300
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogcambia}}

## pEarleyGate 201：另一植物载体的类别对照 {#catalogearley}
```r
sequence <- plasmid_example_pEarleyGate_201
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogearley}}

## pSpCas9(BB)-2A-GFP (PX458)：大型表达载体与文字适配 {#catalogcas}
```r
sequence <- plasmid_example_pSpCas9_BB_2A_GFP_PX458
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogcas}}

## 下一步 {#next}
比较不同密度后，进入[综合实战](plasmids.html#complex)，添加独立引物与酶切图层。
