<a id="anchors"></a>

## 绘制全部特征，只标出选中的记录 {#anchors}

pUC19 示例包含九条参考特征。每条特征有稳定 ID，标签可据此找到对应实体。先用完整表分配道位，再筛出两条 CDS 作标签。实体图层的 ID 是 `features`，标签层通过 `anchor_layer` 指向它。

{{example:feature-subset-labels}}

<a id="placement"></a>

## 决定形状与标签如何放置 {#placement}

混合链向或未知链的区间可用 `geom_feature()`。核酸参考上链向已知的编码特征可选 gene 预设绘制有向箭头。教学表若同时保留原始链和重编码的展示链，应明确传入来源链。若想把多道放在同侧，可给 `position_feature_stack()` 指定 `base_position = position_offset(...)` 和分道方向。稳定特征键与原区间使筛选、重排后的标签仍可跟随道位。

需要指定文字角度时使用固定位置。需要自动排布时使用 `placement = "repel"`，用 `overflow` 决定空间不足时如何处理文字。最后按[目标导出尺寸](output.md)检查标签。
