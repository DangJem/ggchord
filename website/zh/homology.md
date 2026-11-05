<div class="task-intro"><p>比较两条质粒的 AmpR 邻域。先核对真实编码 DNA，再建立对应关系，最后聚焦到 CDS 两侧各 200 bp。</p></div>

## 提取实际 CDS 与编码 DNA {#membership}
选择随包 pUC19 与 pETDuet-1 的 AmpR 注释，保留 ID、坐标、链和文字。按负链进行反向互补，保证比较的是同一生物学方向的 DNA，而非仅比较注释名称。

{{code:homologinput}}

{{output:homologinput}}

把负链 CDS 转成编码方向，再检查两段 DNA 是否完全相同。

{{code:homologdna}}
{{output:homologdna}}

`cluster_features` 应有两条实际 CDS，`cds_dna` 保存各自的定向编码序列。本例的严格目标是完全相同 DNA；一般同源性分析可使用外部比对与聚类结果，不能要求所有同源基因完全一致。

## 显式声明成员和序列对 {#links}
用编码 DNA 的完全一致性建立成员表；不同 DNA 不会因为都叫 AmpR 就分入同一组。

{{code:homologmembership}}

{{output:homologmembership}}

查看 `membership` 与 `homologs$report$comparisons`、`$members`。连接工具消费已有分组，不计算新的序列比对。未知成员与歧义 paralog 会被报告；只有明确选择 `paralogs = "all"` 才请求受界限约束的所有组合。

## 聚焦到可比较的邻域 {#focus}
全长质粒图会把单个 CDS 压得很小。这里从两条真实 CDS 的源范围构建窗口，两侧各增加 200 bp；原始序列表与注释表保持不变。

{{code:homologs}}
{{figure:homologs}}

两条蓝色箭头表示负链 AmpR，灰色带连接已核对的 CDS 区间。刻度来自各自源序列，因此同一个横向位置可以显示不同源坐标。断点标记说明图中只保留了局部窗口。

`local_homologs$report` 记录保留与裁剪情况。窗口外的基因没有显示，不表示它们不存在。区间连接表达分组对应范围，不是逐碱基比对；跨分段注释的外包范围也不等于外显子级对应。

## 将有证据的蛋白区间投影到两端 CDS {#projection}
{{code:linkprojection}}

{{console:linkprojection}}
{{results:projection}}

这里选择相同 CDS 的前 20 个残基，因此 `ungapped = TRUE` 有实际序列依据。普通 BLASTp HSP 不能仅凭两端跨度相等作出这个断言。工具保留源蛋白统计，不把蛋白相似度冒充核苷酸相似度。

## 完成后检查 {#check-result}
两条 AmpR 箭头均沿负链方向。对应关系来自编码 DNA；源坐标与聚焦报告解释了局部显示。保存时使用 [10 × 5 in 的明确尺寸](output.html)。

下一步：[真实 BLASTp 查询区间读取](alignments.html#protein)、[CDS 投影及分段限制](projection.html)、[独立局部窗口](views.html)。

## 保存本页结果 {#save-result}
在工作目录中保存刚刚创建的图；预览时使用相同宽高。

```r
ggsave("ampR-neighbourhood.png", p_homologs, width = 10, height = 5, units = "in", dpi = 144)
```
