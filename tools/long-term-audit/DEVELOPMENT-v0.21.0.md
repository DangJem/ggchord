# v0.21.0 — 真实对齐坐标与局部字符视图

状态：本版约定的本地开发工作完成，不是正式发布验收。2026-11-12用户明确授权将
v0.20.0 b5dccaf快进整合main，并创建v0.21.0持续开发至本版范围完成，已依次执行。
使用原OneDrive checkout，无新worktree。main停在b5dccaf；本版没有整合、push或发布授权。
包开发版本0.21.0.9000；最近正式版仍v0.16.6。

## 实际公开能力

- `as_ggchord_alignment(sequences, alphabet, source_data)`：普通命名字符串或
  accver/sequence表转普通长表，显式DNA/RNA/AA、唯一ID、同宽及合法字母表。
  对齐列保留全部字符/大小写，'-'/'.'消耗列但source_position为NA。
  默认源位置只对应去gap输入片段；不声明完整基因组长度、真实链或同源性。
- source_data显式提供每行的首个非gap start与strand，可指定source_accver。
  负向位置递减，输入字符串须已处于比对方向；可提供完整正向sequence核对，
  核酸负向用互补，AA只反序。source_verified是记录的非gap字符一致性状态，
  在gap行也保留，不表示独立确认来源/基因组注释。
- `read_aligned_fasta()`读取单个已有对齐文件，复用FASTA解析，保留物理头行与可选文件列。
  来源属性区分对齐列和可缺失源位置；不拼接不同对齐，不解析A3M小写插入语义。
- `map_alignment_position()`按对齐行ID在保留单元中精确查找source↔column。
  gap列可匹配但无源位置；窗口外/未知/NA输入为未匹配，不插值、不重编号。
- `geom_residue()`提供独立的标准Cartesian ggplot局部字符层，复用aes/data/position/
  show.legend/inherit.aes、普通样式和用户scale。支持identity/reverse位置标度，
  不接入ggchord生物学骨架Coord或Fusion，不新增Coord。
  draw时按最终panel和字体测量格宽高：可读时画字，过窄按显式tile/error策略处理，
  不缩小字号、不改变源表。图层独立性、数据callback与用户scale优先已直接验证。

运行时首批提交7177c72；后续帮助文字、证据和相关文档收敛另作稳定提交。

## 数据来源与许可

examples/alignment包含Biopython官方Opuntia核酸与蛋白原始CLUSTAL文件、完整许可及来源清单。
仅把每个record的原始blocks拼接成FASTA，ID/大小写/gap/行序不变，不运行任何比对工具。
data-raw/prepare-alignment-inputs.py检查原始SHA256并重建；重建前后两FASTA SHA256一致。
源URL、取得日期2026-11-12、hash、转写规则及Biopython许可均随输入保留。没有复制解析或绘图实现。

DNA有7记录×906列=6342单元，蛋白20记录×411列=8220单元，共1295个gap单元。
原始来源位置默认是去gap片段的坐标，不是完整基因组坐标；负向验证是对原对齐明确重编码的
方向控制，不冒充新的分子链证据。source_data演示从同一对齐去gap重建的字符串只验证内部
一致性，教程明确不当作独立基因组验证。

## 定向验证

`Rscript tools/long-term-audit/check-alignment-detail.R`通过：

- 两组真实输入所有记录的非gap源位置与原始字符往返、gap保留、文件/物理头行来源；
- 明确声明的负向重编码、IUPAC互补、完整正向字符串核对，以及已知来源不一致拒绝；
- 标明合成的RNA/大小写与基础格式错误，窗口/未知/NA查找和gap列的匹配区别；
- 三条DNA的原始321–380列与三条蛋白221–260列，保留真实gap和字符差异；
- 实际draw时宽画布显示字母、窄画布明确报错或保留色格；reverse位置标度通过，
  非Cartesian/非线性位置标度拒绝；两个独立图层和显式用户scale通过；
- 源数据不变、临时设备清理恢复。R在预期的draw-time错误设备关闭时可能报告
  `Killing locked device`；记录了这个清理warning，设备列表恢复断言仍通过。
  普通成功PNG无此warning，不把可选错误探针的warning当成字符被成功绘出。

输出v0.21.0-alignment-results.json、两组窗口CSV与DNA/蛋白/窄格三张PNG。
最终两张字符PNG和窄格图已人工查看；本机耗时仅为该脚本记录，不设性能门槛。
原import/focus测试不扩张；只添加一个覆盖适配、gap查找和基础gtable的最小成功路径，未运行testthat。

## 文档与验证范围

NEWS、四个公开帮助、汉语参考、数据契约、两语formats课程、兼容矩阵、安装版本、
API/参数/任务索引同步；网站输入包生成规则将新输入和许可纳入未来构建。
145导出、1592形参、109 Rd、1022段参考翻译。formats在独立R环境执行准备及raw R代码通过；
双语代码/锚点一致，124源页面的上下文lint、教学索引与签名清单通过。
已同步Rmd源，但没有编译vignette。

**整站 source/assets 总检查没有通过**：5个旧lesson engine hash和178张旧图的保守依赖指纹已落后。
这些涉及v0.20.0以来的源变更，不能用新字符脚本替代旧教程执行或把旧图stamp成当前输出。
明细在v0.21.0-documentation-scope.json；这是整站重新验证待办，不宣称整站/浏览器验收完成。
本轮只核对本版新图、相关示例和源文档，不重建冻结docs、.site或全套旧PNG。

未运行完整testthat、devtools::test、R CMD build/check、vignette编译、PDF/SVG/PNG发布矩阵。
README未来方向、原有man/figures、原始旧输入与许可证保持不变。第三方新输入自带许可。

## 未纳入本版的后续能力

CIGAR/cs/Clustal/Stockholm/A3M全格式解析、MSA计算、共识/Logo、自动CDS翻译/phase、
宏观字符投影、公开magnify API和M1–M3、Fusion F2/F3、通用Guide/Adaptive及1.0治理
仍按v10各自门槛推进。本版完成的是已有对齐→精确源映射→可编辑局部字符图这一可用工作流，
不把它扩大为整个MSA或Semantic Enrichment完成。
