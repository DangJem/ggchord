<div class="home-intro">
<img class="package-logo" src="../assets/logo.webp" alt="ggchord" width="320" height="320">
<p class="eyebrow">GENOMIC VISUALIZATION · R / GGPLOT2</p>
<p class="home-lead">用 R 绘制基因组比较和质粒注释</p>
<p>用图层组合弦图、共线性图与密集质粒注释。从真实序列、比对和注释，构建可复现的科研图形。</p>
<div class="home-actions"><a class="primary-action" href="installation.html">开始使用 →</a><a href="about.html">浏览章节</a></div>
</div>

## 选择你的可视化任务 {#choose}
- [弦图](chord-start.html) · 比较多条序列的对应区间。
- [共线性图](collinear-start.html) · 按共享比例比较基因位置与顺序。
- [质粒图](plasmid-start.html) · 在单条序列上标记功能区间。

## 弦图 · 四条基因组的总体联系 {#chord}

{{figure:showcasechord}}

红色正链箭头放在内侧，连接带在重叠的端点局部让位；蓝色负链箭头放在外侧。连接带颜色表示比对相似度。两条 MT 序列的原始注释只有正链 CDS，另外两条同时展示正负链。[数据、筛选与完整代码](chord.html)。

## 共线性图 · 四序列、多图层、跨行比较 {#collinear}

{{figure:showcasecollinear}}

沿同一碱基尺度比较基因顺序与对应位置。上下箭头区分链方向，虚线轮廓标记跨行连接。[共线性图入门](collinear-start.html) · [完整案例](collinear.html)。

## 质粒图 · pETDuet-1 的密集注释 {#circular}
{{figure:showcaseplasmid}}

双链骨架、多层特征、五个参考引物，以及按同一切点合并的酶名称。图中输入包含 58 个切点和 20 个特征，文字按实际导出尺寸布局。
[从简单圆图到完整质粒注释 →](plasmid-start.html) · [完整案例](plasmids.html#complex)

## 按你的起点继续 {#learn}
- [导入你自己的数据](formats.html)
- [调整配色、图例和主题](styling.html)
- [查看全部教程与功能](about.html)
- [查阅函数、参数和返回值](reference/index.html)
