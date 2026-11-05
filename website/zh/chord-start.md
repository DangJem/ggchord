从随包的四条序列开始。先看各条骨架的位置和长度，再添加比对和特征。

## 绘制骨架 {#draw}

{{code:frontbackbone}}

{{figure:frontbackbone}}

## 读图 {#read}

每条弧代表序列表的一行。用标签找到 `MT108731.1`，再与 `MT118296.1` 比较。它们分别为 64,323 bp 和 32,090 bp，所以第一条弧约为第二条的两倍长。

`geom_seq()` 画出骨架，`geom_seq_label()` 标出序列 ID，`coord_chord()` 把四条骨架排成弧形。

## 试着去掉标签 {#try}

删除 `geom_seq_label()` 那一行并重新运行。四条弧会保留，序列名会消失。

## 换成自己的序列 {#adapt}

把 `seq_data_example` 换成含 `accver` 和 `length` 的表即可。列的含义见[基本数据表](chord-data.html#tables)。接下来[添加连接和特征](chord-build.html)。
