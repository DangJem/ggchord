## 范围与时间 {#scope}
下列写法是当前候选替代，不是 v1.0.0 长期兼容承诺。0.16.0–0.16.4 阶段成果已随 v0.16.5 交付。
当前 `0.16.6.9000` checkout 尚未发布，正式安装标签仍为 v0.16.5。RC 审计保留
灵活参数、拼写别名、重导出及 getter/export 入口，没有新增弃用。
表中的旧调用用于解释迁移，不应作为当前例子执行。
历史证据只提供删除上界时，不猜测精确首删 tag。

## 替代表 {#table}
{{migration}}

## 运行当前替代 {#replacement}
在新会话中运行下面的准备代码，再执行当前替代写法：

{{prepare:setup,plasmiddata}}

{{code:migration}}

旧负链 offset 幅度在显式 Position 中需要负号。映射仍使用
`aes(feature_shape = type)` 配合 shape 标度，固定参数则为 `shape`。
自动标签使用 `label_layout`，手动几何控制使用相应固定标签层。
`label_layout = "aligned"` 改为 `"auto"`。已删除名称明确报错，不保留静默执行路径。

## 保留的支持 {#retained}
五个 colour/color 拼写别名、全部上游重导出、gene/注释/plasmid 快捷入口及 getter/export
两类几何接口均保留。各坐标默认差异和已记录结果字段保持。NULL/NA/缺参语义应查阅
[接口说明](reference/ggchord_api_contract.html)和当前参考；NULL 不统一表示关闭。
