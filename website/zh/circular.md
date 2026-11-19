单序列圆图把一条参考序列的位置排在圆环上。输入为序列全长和注释表。各特征与骨架使用同一坐标单位，序列首尾在圆环接缝相遇。

<a id="draw"></a>

## 画出 pUC19 的九条注释 {#draw}

随包 pUC19 记录提供 2,686 bp 的完整长度和九条已有注释。`geom_feature()` 直接读取特征表；已知方向的箭头沿注释链绘制。

{{example:first-circular}}

<a id="adapt"></a>

## 加入标签与位点 {#adapt}

给 `geom_feature()` 和 `geom_feature_label()` 传入同一位置调整，可使标签跟随特征所在道。引物区间和酶切位点各用独立图层。完整 pUC19 记录含 39 条位点记录，对应 33 个不同位置；同一位置可能有多个酶名。

[pUC19 完整注释图](case-circular.md) · [十三张质粒图](plasmid-catalog.md) · [圆形坐标](../../man/coord_circular.Rd)
