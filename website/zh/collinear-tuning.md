用原来的四条骨架和一条比对作对照，分别调整行距和加入标尺。

## 默认布局 {#baseline}

{{code:frontlinearbase,frontlinks,frontlinearlink}}

{{figure:frontlinearlink}}

## 增加行距 {#gap}

{{code:frontlineargap}}

{{figure:frontlineargap}}

`seq_gap` 在共线性图中以最长骨架宽度为基准，默认是 `0.30`。只把 MT108731.1 下方的间隙设为 `0.45`，其余行距仍为 `0.30`。比对带仍连接相同的 bp 区间。

## 添加固定长度标尺 {#scale}

{{code:frontlinearscale}}

{{figure:frontlinearscale}}

回到默认行距，加上 10 kb 标尺。用它对照短序列的长度；它约占 32 kb 骨架的三分之一。下一节介绍[输入和显示方向](collinear-data.html)。
