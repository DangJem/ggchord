这些调整属于不同职责：geom控制逐序列几何，Position控制局部放置，coord控制整体视图，scale与guide负责颜色和图例。

## 本页准备 {#prepare}
先解压[教学输入](../assets/tutorial-inputs.zip)。
{{prepare:sequencepalette,docstyle,showcasedata}}

## 先检查方向证据 {#orientation}
`show_optimization` 只优化显示方向，保留序列顺序；按真实比对长度加权。检查其seq_orientation与score_before/score_after，不能把显示翻转解释为新的生物学倒位结论。自动优化是启发式，不保证全局最优。

## 将功能组合成可读的图 {#parameters}
{{code:showcasechord}}
{{figure:showcasechord}}
命名向量逐序列匹配：曲率0是直线，正负值向相反方向弯曲，绝对值影响弓高；半径改变序列位置。上图同时演示自动方向、不同半径、零/正/负曲率和局部避障。参数不修改源坐标。

## 避障方式为何看起来不同 {#clearance}
`link_avoid = "uniform"` 在整个端区间留统一空隙；`"smooth"` 只在特征附近平滑退让，形成局部缺口；`"none"` 不做特征避障。`link_gap = NULL` 才启用自动规则，给出固定数值会覆盖自动间距。避障检测箭头/特征几何，不避让文字，也不会自动绕过跨行连接途经的整条中间序列。

## 颜色与导出 {#guides}
identity采用独立左侧长色条，序列和链图例在底部。链图例显式传递方向覆盖，箭头符号与图中几何一致。[样式专题](styling.html)解释图例职责；按12×10 in导出，最后到[1.7 首页综合示例](chord.html)核对数据与结果。
