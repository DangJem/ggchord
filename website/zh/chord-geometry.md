{{code:courseblockb}}

{{prepare:frontbackbone}}

## 旋转整张图 {#orientation}

{{code:frontrotation}}

{{figure:frontrotation}}

`rotation` 控制逆时针旋转角度，默认是 0°。改为 `rotation = 90` 后，整组骨架和标签相对默认图逆时针转过 90°。对照默认图，观察相邻的弧之间的空隙移到了哪里。旋转可以给某一侧腾出排版空间，序列内部的坐标顺序保持原样。

## 加大骨架间隙 {#parameters}

{{code:frontgap}}

{{figure:frontgap}}

`seq_gap` 控制相邻骨架间的空白比例。默认是 `0.03`，本例只把 MT108731.1 之后的间隙提高到 `0.12`，其余仍用 `0.03`。 最后一条序列之后的间隙位于环的首尾接缝。和默认图对照时，先看相邻的弧之间的空白，再看相邻的弧仍保有长短之别。

## 半径：向内移动一条弧 {#shape}

{{code:frontradius}}

{{figure:frontradius}}

默认半径为 `1`。另外三条弧保留默认值，`MT118296.1` 用 `0.72`。对照[默认骨架](multi-sequence-start.html#draw)，看这条弧怎样向图心收进，序列末端仍在约 32,090 bp。半径是布局单位，读序列长度要看刻度。

## 曲率：把弧画成直线 {#curvature}

{{code:frontcurvature}}

{{figure:frontcurvature}}

`seq_curvature = 1` 是标准圆弧，`0` 把该骨架画成直线，负值会向相反方向弯。其余三条骨架保留标准圆弧。若同时调整半径和曲率，先分别观察两个参数，再把它们合用。
