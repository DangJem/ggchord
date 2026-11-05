用相同的四条序列，先加一条比对，再加一个 CDS。本页代码依次运行即可。

## 重建四条骨架 {#bundled}

<span id="prepare"></span><span id="switch"></span><span id="backbones"></span>

{{code:frontlinearbase}}
{{figure:frontlinearbase}}

## 加入一条区间比对 {#links}

{{code:frontlinks,frontlinearlink}}

{{figure:frontlinearlink}}

连接带在第一行覆盖 26,298–51,139 bp，在第二行覆盖 7,121–31,959 bp。横向错开体现两个区间在各自序列上的位置。

## 查看 CDS 并加上箭头 {#inspect-genes}

{{code:frontgenes,frontselectedgene}}

{{output:frontselectedgene}}

这条原始正链 integrase CDS 位于第三行的 43,765–46,140 bp。表中 `anno` 保留原始注释名。

### 加到图上 {#add-genes}

{{code:frontlinearcds}}

{{figure:frontlinearcds}}

<span id="labels"></span>红色箭头位于第三条绿色骨架的右半部。需要显示名称时，可另加 `geom_gene_label(data = chord_genes_example)`。接下来[调整行距与标尺](collinear-tuning.html)。
