用相同的四条序列，每次只改一个参数。先看默认图，再比较骨架的位置与留白。

## 默认图 {#prepare}

{{code:frontbackbone}}

{{figure:frontbackbone}}

## 旋转整张图 {#orientation}

{{code:frontrotation}}

{{figure:frontrotation}}

`rotation` 控制逆时针旋转角度，默认是 0°。改为 `rotation = 90` 后，整组骨架和标签相对默认图逆时针转过 90°。对照默认图，观察相邻的弧之间的空隙移到了哪里。旋转可以给某一侧腾出排版空间，序列内部的坐标顺序保持原样。

## 加大骨架间隙 {#parameters}

{{code:frontgap}}

{{figure:frontgap}}

`seq_gap` 控制相邻骨架间的空白比例。默认是 `0.03`，本例只把 MT108731.1 之后的间隙提高到 `0.12`，其余仍用 `0.03`。 最后一条序列之后的间隙位于环的首尾接缝。和默认图对照时，先看相邻的弧之间的空白，再看相邻的弧仍保有长短之别。

## 什么时候需要更多调整？ {#clearance}

比较图里如果出现连接带和箭头重叠，查看[连接避障](chord-routing.html#clearance)；如果需要给某条骨架单独调整半径或弧度，查看[逐序列参数](parameters.html)。

<span id="guides"></span>下一节[逐条调整半径、曲率和显示方向](chord-geometry.html)。
