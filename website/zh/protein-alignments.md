使用两份内置表比较四条完整蛋白，无需先下载文件或运行 BLASTp。蛋白骨架与连接端点均使用氨基酸坐标（aa）。

## 直接绘制 BLASTp 蛋白比对 {#protein}

`protein_seq_data_example` 包含三条 EMALE 蛋白和一条 Mavirus 蛋白。先查看标识、完整长度和单位：

{{code:nativeproteininput}}
{{output:nativeproteininput}}

长度分别为 302、302、306 和 310 aa，来自完整蛋白序列，不是命中的最大端点。`sequence` 列还保留完整氨基酸序列，这里只展示绘图需要的三列。

`protein_ribbon_data_example` 包含四条蛋白之间六对实际 BLASTp 比对。下面查看两端范围及一致性：

{{code:nativeproteinlinks}}
{{output:nativeproteinlinks}}

{{code:nativeproteinplot}}
{{figure:nativeprotein}}

共线性布局默认连接相邻行，显示三条带；颜色表示一致性，下方标尺长 50 aa。前两条蛋白高度相似，后两对相邻比较更弱。表中区间包含两端，含缺口的比对保留两个端点范围。

改为弦图布局时，默认显示全部六对连接，可同时比较任意两条蛋白：

{{code:nativeproteinchord}}
{{figure:nativeproteinchord}}

六对比对的一致性约为 35.7%–99.0%。图中颜色沿用默认的一致性标度，无需另外设置配色。

蛋白区间注释使用 `geom_feature()` 等通用图层，以氨基酸位点提供范围，并用 `directionality = "none"` 声明无方向区间。同一张图统一使用 aa 或 bp；`blastx`、`tblastn` 的两端单位不同，应先转换到统一坐标，或分别展示。基因箭头、引物与酶切图层用于核酸输入。圆形蛋白布局仅表示排列方式。

## 将蛋白命中显示在基因组上 {#projection}

如果需要在核酸骨架上定位命中，则提供已知 CDS 映射，将 aa 区间投影到 bp。下面另用 BVI_023A_0009 的源蛋白命中演示，先选出具有完整 CDS 的记录，再转换 query 端。此小节读取[说明输入包](formats.html)中的原始文件。

```r
library(dplyr)
```

{{code:protein}}
{{output:protein}}
{{code:proteinproject}}
{{output:proteinproject}}
{{results:protein}}

氨基酸范围 50–289 投影后对应 BVI_023A 的 7,181–7,900 bp，方向为负链。单侧范围投影可保留有缺口 HSP 的覆盖范围。需要连接两端基因组时，还应提供 subject 的 CDS 模型；双端无缺口投影见[CDS 投影](projection.html)。

## 数据来源 {#source}

内置蛋白标识为 BVI_029_0020、BVI_058_0011、E4-10_023_0015 和 ADZ16414.1。完整 FASTA 来自公开 EMALE/Mavirus 数据，采用 NCBI BLASTp 2.17.0+ 实际计算全对全比对，每对保留一个 HSP；原始文件保留自身命中及两个比较方向，共 16 行。说明表按序列表顺序为六对各保留一个方向，未修改端点或一致性。

输入包的 `examples/protein/` 保留完整 FASTA、原始 BLASTp 结果及 `sources.json` 来源哈希；数据生成规则见 `data-raw/prepare-bundled-proteins.py` 和 `data-raw/generate_protein_data.R`。绘图直接使用内置对象，无需安装 BLASTp。

## 检查要点 {#check-result}

核对完整长度、命中端点与单位。原生蛋白图以 aa 为轴；定位到基因组时使用 CDS 模型投影到 bp。保存前用 `view_ggchord()` 检查目标尺寸。
