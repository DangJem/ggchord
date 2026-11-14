<div class="task-intro"><p>把已知蛋白区间投影到真实 CDS 的基因组位置，显式保留 phase、链与分段关系。</p></div>

## 本页准备 {#prepare}
准备[已核对编码 DNA 的 AmpR 数据](homology.html)，或将下方脚本放入工作目录并运行。

{{recipe:ampr}}

## 先准备有来源的 CDS {#cds}
本页准备代码包含[完全相同 CDS 的比较](homology.html#membership)，使用随包质粒的真实注释
和 DNA。完整 CDS 的 phase 明确为 0，翻译表显式设为 11。

{{code:projection}}

{{console:projection}}

这个例子选择第 2–20 个残基作显示范围，不声称它是检测到的蛋白结构域。
投影结果是基因组区间，保留正负链与分段关系。真实检测到的 domain 也应提供确切
父 CDS ID 及氨基酸范围；不能通过同名注释或注释最大端点猜测。

## 分段、phase 和部分区域 {#boundaries}
`build_ggchord_cds()` 先建立 CDS 分段映射，`project_ggchord_features()` 再将氨基酸
每个输入 feature 返回一行，跨剪接片段放在 `segments` 列表列；
`project_ggchord_links()` 才按双方边界拆成配对行。不要简单把坐标乘以 3。
不完整 CDS、phase、截断或不支持的位置必须按文档策略处理并检查报告。
DNA 匹配、氨基酸区间投影和真实翻译/比对是不同步骤，绘图库不会代为推断缺失证据。

## BLASTp 查询范围与双端连接 {#links}
{{code:linkprojection}}

{{console:linkprojection}}

这个无缺口例子依赖[同源专题](homology.html)中已核对的相同编码序列。
真实 BLASTp 输入可以有 gap；查询区间仍能投影，但双端连接还需要 subject 的 CDS
和无缺口对应证据。请阅读[真实 BLASTp 读取与投影](protein-alignments.html#projection)。
绘制以 bp 为轴的图时，应使用投影后的基因组端点。

## 检查要点 {#check-result}
检查父 CDS ID、投影报告及返回分段。选择的残基区间只是显示例子，不能将其描述为检测到的结构域。

<div class="next-step">继续：<a href="alignments.html">导入核酸和蛋白比对 →</a></div>
