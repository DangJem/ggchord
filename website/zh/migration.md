## 范围与时间 {#scope}
下列写法是当前候选替代，不是 v1.0.0 长期兼容承诺。此前 0.16.x 版本均本地完成、
尚未发布。表中的旧调用用于解释迁移，不应作为当前例子执行。
历史证据只提供删除上界时，不猜测精确首删 tag。

## 替代表 {#table}
{{migration}}

## 运行当前替代 {#replacement}
先运行[准备代码](index.html#quick-start)及[圆图](modes.html#circular)中的 `features` 表：

{{code:migration}}

旧负链 offset 幅度在显式 Position 中需要负号。映射仍使用
`aes(feature_shape = type)` 配合 shape scale，固定参数则为 `shape`。
自动标签使用 `label_layout`，手动几何控制使用相应固定标签层。
`label_layout = "aligned"` 改为 `"auto"`。已删除名称明确报错，不保留静默执行路径。

## 保留的支持 {#retained}
五个 colour/color 拼写别名、全部上游重导出、gene/feature/plasmid 快捷入口及 getter/export
两类几何接口均保留。各坐标默认差异和已记录结果字段保持。NULL/NA/缺参语义应查阅
[API 契约](reference/ggchord_api_contract.html)和当前参考；NULL 不统一表示关闭。
