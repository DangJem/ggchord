{{sitehero}}

## 选择你的可视化任务 {#choose}

{{taskcards}}

## 从零开始，只需三步 {#quick-start}

1. 按[安装页](installation.html)安装包并确认所需版本。
2. 加载 ggchord，先使用随包真实序列数据。
3. 选择布局；有图之后再加连接、注释和标签。

```r
library(ggchord)

p <- ggchord(seq_data = seq_data_example) +
  geom_seq() +
  geom_seq_label() +
  coord_chord()

view_ggchord(p, width = 9, height = 6, units = "in")
```

同一套序列表可以采用不同布局；单序列圆图需要只显示一条序列。各类[独立首课](learn.html#modules)分别解释输入、读图与下一步。

## 真实数据，完整结果 {#results}

以下图形来自当前 ggchord 源码和可追溯输入。预览用于选择案例；点击可查看高清原图。

<div class="home-results">

<section>

### 四基因组弦图 {#chord}

{{figure:showcasechord}}

[来源、选择与完整代码](chord.html)

</section>

<section>

### 四基因组共线性图 {#collinear}

{{figure:showcasecollinear}}

[来源、选择与完整代码](collinear.html)

</section>

<section>

### 密集质粒注释 {#circular}

{{figure:showcaseplasmid}}

[来源、选择与完整代码](plasmids.html#complex)

</section>

</div>

## 按你的起点继续 {#learn}

- [导入自己的数据](formats.html)
- [按任务找函数与参数](learn.html#tasks)
- [调整颜色、图例与标签](styling.html)
- [按论文尺寸预览与导出](output.html)
- [查阅完整 API](reference/index.html)

## 复现你的图 {#versions}

导出时保存包版本、输入文件和绘图代码。完整案例提供数据来源与整理取舍；版本要求见[安装页](installation.html)。
