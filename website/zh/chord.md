想了解多条序列中哪些区间互相匹配时，可以画弦图。弧线表示序列骨架，带子的两端各占据来源序列上的一段区间。

<a id="draw"></a>

## 在四条骨架上画一条比对 {#draw}

`ribbon_data_example` 首行连接四条参考序列中的两条。query 在 MT108731.1 上覆盖 26,298–51,139 bp，subject 在 MT118296.1 上覆盖 7,121–31,959 bp。另外两条弧保留在图中，便于比较全长。

{{example:first-chord}}

<a id="adapt"></a>

## 换用自己的比对 {#adapt}

每条骨架提供 `accver` 和完整 `length`。连接表需要 `qaccver`、`saccver` 以及 `qstart`、`qend`、`sstart`、`send` 四个端点。反向比对产生的降序端点应保留。要用色标比较多条带的一致性，还需保留连接表中的 `pident`。

[四基因组弦图](case-chord.md) · [区间连接层](../../man/geom_link_ribbon.Rd) · [弦图坐标](../../man/coord_chord.Rd)
