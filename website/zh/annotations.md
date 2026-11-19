当圆图上的注释较密时，通常希望保留全部特征，只给其中一部分加标签。可以先用完整的特征表安排布局，再按记录中的字段选择要标注的条目。

<a id="placement"></a>

## 先确定特征位置，再选择标签 {#placement}

在 2,686 bp 的 pUC19 骨架上画出九条注释。每一行特征都有登录号和起止位置，方向已知时也记录链别。靠近序列起点的一组短区间有多处重叠。

{{example:first-circular}}

<a id="anchors"></a>

## 保留九条特征，只标注两个 CDS {#anchors}

保留源记录的匹配标识作为 `feature_id`。将实体层命名为 `features`，用 `position_feature_stack()` 把重叠区间分到不同的道上。标签层只取两个 CDS，再通过 `anchor_layer` 指定原实体层。这样筛选之后，每个标签仍跟随各自的特征。

{{example:feature-subset-labels}}

在起点附近可以找到 lacZα，圆图另一侧较长的区间标有 AmpR。其余七条记录仍然可见。下一步可[按注释类型设置配色](styling.md)，也可以查看[完整 pUC19 图的位点标记](case-circular.md)。
