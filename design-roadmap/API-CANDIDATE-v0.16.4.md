# v0.16.4 API 候选冻结与迁移基线

**状态：本地开发完成；尚未正式验收或发布。**

用户于2026-10-31明确启用本版。冻结对象为当前公开合同，供v0.16.5文档和v0.16.6
候选审计共同引用；不等于已发布或v1.0.0长期兼容承诺，不冻结私有算法和浮点布局。

## 唯一当前清单

| 内容 | 当前事实与来源 |
| --- | --- |
| 公开入口 | `candidate-exports-v0.16.4.csv`：130个导出，108个自有函数、22个原样重导出；5个直接拼写别名保留 |
| 参数与默认 | `candidate-parameters-v0.16.4.csv`：1,236个自有参数，当前formals/Rd、单位、NULL/NA/缺参、职责和旧语义来源分列 |
| 返回 | `candidate-returns-v0.16.4.csv`：130项当前返回说明；数据34入口和输出16组字段继续引用其详细合同 |
| aesthetic | `candidate-aesthetics-v0.16.4.csv`：原绘图角色及当前feature/track/通用样式/拼写角色；明确输入/后期阶段和适用入口，不冻结私有draw字段或scoped key编码 |
| 坐标支持 | `candidate-compatibility-v0.16.4.csv`：继承69项源码支持合同及有限运行出处，不把NOT_RUN改成PASS |
| S3 | `candidate-s3-v0.16.4.csv`：9项注册保留；上游派发协议不等于内部ggproto字段稳定 |
| 迁移 | `candidate-migrations-v0.16.4.csv`：200条参数、取值、映射/字段及旧函数记录；当前替代、删除上界、首删版本证据与历史记录分开 |
| 来源 | `candidate-source-hashes-v0.16.4.csv`：源码、Rd、旧审阅输入与生成器MD5；HEAD和实际工作树字节明确区分 |

自有参数语义复用已有1019项绘图/基础审阅、195项v0.16.2数据和22项v0.16.3输出
审阅。重命名和七个共享标签参数按实际共享内核定位，当前Rd覆盖历史名字的文案；
保留旧函数/参数名与审阅SHA，不将全部参数升级为本版分别运行通过。旧shape简写
曾允许NULL并让映射优先，该语义已于v0.16.1被唯一固定shape入口取代；本版整合表
显式覆盖该历史说明和ribbon alpha的旧优先描述，缺参说明按当前默认重建。

## 已收束的公开取舍

- F07：getter/export均保留已描述成员。export是分层表格推荐入口，getter提供完整
  布局；preview返回临时图片路径，diagnose返回目标尺寸的有限测量证据。
- F10–F13：固定shape、width、label_*、link_gap/link_ctrl_point为唯一当前参数名；
  mapped feature_shape保留。旧offset于v0.16.1删除，显式逐链Position值带符号。
- F14/F15：5个colour/color直接别名与22个ggplot2/grid重导出全部保留；文档首选
  colour。重导出保持上游对象身份，其未来上游签名/内部结构不由ggchord冻结。
- F24/F28：保留已删除入口的迁移错误，不恢复兼容执行。constructor拒绝未知dots；
  普通geom/scale/guide/theme仍按各自已记录的上游委派规则处理额外参数。
- F51：自动标签层的六个手动label_*控制原先被忽略，现明确指向对应固定标签层。
  未知gene_label_*不再建议不存在的新参数，aligned直接指向label_layout="auto"。
  没有删除当前公开形参或改变有效默认/标准文本样式。

功能职责和默认差异保留：gene是严格定向快捷入口，feature容纳未知链；plasmid预设
有独立默认。line默认avoid=none，ribbon默认uniform；chord的NULL控制点仍为原点。
三种坐标的支持边界不统一成假兼容，基因组锚点历史显示公式不在本版改动。
固定样式优先，省略固定值保留映射，用户scale优先，同类层各自拥有data/mapping/参数。

## 候选合同的含义

formals中的NULL是语法默认，实际fallback继续由对应入口、theme/preset和coord解释；
不能把所有NULL都解释成关闭或缺失。bp、layout data、角度、数值域和物理mm分别定义。
clean sort/prepare keep、file/files、未知链、单bp方向、来源行和聚合估计标记保持。
数据转换/辅助搜索的能力边界见 `ggchord_data_contract`；结果字段与缓存生命周期见
`ggchord_layout_contract`，跨域入口说明见新增帮助 `ggchord_api_contract`。

已在Rd承诺的字段保持：annotation_registry、sequence_layout、conversion、track、
label/connection诊断等不能通过改称“内部”删除。source_row是resolved input行，
不是统一物理文件行；装饰可为NA，聚合另有source_rows。资源号不承诺跨编辑永久稳定。
测量panel与geometry收敛分开，字号和布局可以随设备/字体改变。

## 迁移的当前执行方向

| 旧调用 | 当前调用与影响 |
| --- | --- |
| `geom_feature(feature_shape="block", feature_width=.05)` | `geom_feature(shape="block", width=.05)`；映射仍用aes(feature_shape=...) |
| `geom_gene(gene_offset=c("+"=.1,"-"=.2))` | `geom_gene(position=position_strand(c("+"=.1,"-"=-.2)))`；负链旧幅度须取负 |
| `geom_gene_label_repel(gene_label_layout="aligned")` | `geom_gene_label_repel(label_layout="auto")`，一步完成名字与取值迁移 |
| 自动层的`label_rotation=90` | `geom_gene_label(label_rotation=90)`或`geom_feature_label(...)`；需要自动布局时删除手动控制并选择label_layout |
| `leader="trunk"` | `leader="radial"`；radial默认不改 |
| `ggchord::geom_ribbon()` | `geom_link_ribbon()`；不是删除ggplot2的同名geom |

前述旧语法用于解释迁移，不作为当前可执行示例。当前完整教学输入与替代组合见
`tools/validate-api-candidate.R`和`help("ggchord_api_contract")`。历史API-MIGRATION
小节按各版当时事实保留；特别是v0.16.0的feature_shape示例已被v0.16.1取代。

## 冻结检查与修改政策

`Rscript tools/freeze-api.R check`核对自有公开签名/默认、拼写绑定、重导出原样身份
和S3注册，不绘图，不运行testthat。上游重导出的签名仅记录快照，不把依赖升级误报为
本包形参破坏。该检查不是语义/视觉认证，也不锁定私有字段、文字框或迭代数。

明确决定修改API后，先记录事实、冲突/理由、旧新调用和删除/变更版本，再显式运行
`Rscript tools/catalog-api-migrations.R`与`Rscript tools/freeze-api.R snapshot`更新候选。
默认check不覆盖快照。语义变化、有效默认或返回字段变化仍须人工审阅；源码hash
只是来源指纹，不是行为证明。未确认的想法不得混入候选承诺。

## 后续边界

v0.16.5使用此清单重构英语/汉语教程、README、网站与图片生成源，修复历史教程
过期调用；它可改善文档，不因文案修改重开已完成的API取舍。网站部署与资产政策
按该版完整文档和用户授权执行，本版不重建冻结docs、不提前创建下一版分支。

v0.16.6按实际风险处理低频等价格式、复杂几何、设备/字体/平台兼容和候选发布阻断。
历史删除的“removed_by”与“first_removed”证据区别须在迁移文档保留；没有首删tag
证明时不凭类别NEWS猜精确首删日期。当前不支持/未知范围有明确说明，不冒充通过。
正式发布仍需另行明确授权；本轮不执行testthat、包构建/check、vignette或设备发布矩阵。
