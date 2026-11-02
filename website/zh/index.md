<div class="home-intro">
<img class="package-logo" src="../assets/logo.webp" alt="ggchord" width="320" height="320">
<p class="eyebrow">GENOMIC VISUALIZATION · R / GGPLOT2</p>
<p class="home-lead">结构、关系与注释，一图呈现。</p>
<p>用图层组合弦图、共线性图与密集质粒注释。从真实序列、比对和注释，构建可复现的科研图形。</p>
<div class="home-actions"><a class="primary-action" href="installation.html">绘制第一张图 →</a><a href="articles.html">按任务查找教程</a></div>
</div>

## 选择你的可视化任务 {#choose}
{{showcases}}

## 弦图 · 四条基因组的总体联系 {#chord}
{{figure:showcasechord}}

从ggchord最早的弦图工作流出发：直线与正负曲率骨架、不同半径、自动显示方向、局部ribbon避障与逐一标注的精选特征。比对来自真实观测。
[数据、完整代码与读图说明 →](chord.html)

## 共线性图 · 四序列、多图层、跨行比较 {#collinear}
{{figure:showcasecollinear}}

弦图的直线坐标表达。复用真实序列、比对和相似度配色，上下分列正负链特征，局部避障并向右绕行跨序列ribbon。共同图层可替换坐标；当前圆形专属参数的限制在教程中说明。
[完整流程与各层含义 →](collinear.html)

## 质粒图 · pETDuet-1 的密集注释 {#circular}
{{figure:showcaseplasmid}}

双链骨架、多层特征、五个参考引物，以及按同一切点合并的酶名称。图中输入包含 58 个切点和 20 个特征，文字按实际导出尺寸布局。
[从简单圆图到完整质粒注释 →](plasmids.html#complex)

## 按你的起点继续 {#learn}
- [安装并画第一张图](installation.html)
- [导入你自己的数据](formats.html)
- [调整配色、图例和主题](styling.html)
- [查看全部教程与功能](articles.html)
- [查阅函数、参数和返回值](reference/index.html)

示例优先使用 dplyr、purrr、tidyr 和原生 `|>`；图层组合遵循 ggplot2。需要更新旧代码时查看[迁移指南](migration.html)。
