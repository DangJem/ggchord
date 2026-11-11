# v0.20.0 Semantic Enrichment 开发记录

状态：本版约定的本地开发范围已完成：BEDPE、显式FASTA字符输入、聚焦修复、字符/M0隔离原型和适配器审计。
这不表示整个跨版本Semantic Enrichment阶段或v10长期任务书完成；没有公开字符/语义放大API。
2026-11-11用户明确要求将v0.19.0并入main并开展下一阶段。已将a961b0e快进整合本地main，
再从整合后的main创建v0.20.0；沿用原OneDrive checkout，无worktree。
包开发版本0.20.0.9000；最近正式版仍为v0.16.6。不push、部署或正式发布；本版整合main仍须另获授权。

## 首个实现切片

- read_bedpe()读取6个必填字段、顺序可选的4个标准字段与任意扩展列；每文件字段数一致，严格按tab分隔。
- 第一锚点为query，第二为subject，不排序或交换角色；零基半开坐标只转换一次。
  负链按既有配对区间适配器采用降序canonical端点；未知链不变成正链或已知相对方向。
- 原始标准字段全部保留为字符，score不强制数值化；额外列在bedpe_extra中按物理列号保存，包含尾部空字段。
  记录物理.source_line、可选.source_file、输入/输出坐标与逐文件行数；name不冒充唯一HSP身份。
- 默认拒绝.染色体/-1坐标等不完整锚点；显式unmapped="keep"保留用于审计的缺失canonical字段与q_mapped/s_mapped。
  保留不等于可绘制，绘图前须筛选完整锚点；非法已知坐标和零宽始终拒绝。
- 不推断pident、比对长度、支持度、完整序列长或逐字符对应；不新增专用Geom或外部分析工具依赖。

规范来源：https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format
只参考公开格式描述，解析实现独立编写，没有复制第三方代码或添加外部资产。

## 定向验证

`Rscript tools/long-term-audit/check-bedpe.R`通过：
仓库真实EMALE PAF 11行的明确BEDPE转写，q/s四端点往返一致、额外证据/原始字段/物理行号保留；
一个使用固定填充的Collinear完整Ribbon构建成功，不制造Identity。
转写中的query '+'表示选择的比对参考轴，不是基因链注释或真实配对测序证据。
合成格式用例覆盖六列单位置、未知链、不完整锚点审计、任意score、尾部空额外字段、
多文件不同列数、保留track/browser作为合法tab分隔accession、零宽/坏字段数拒绝。
结果见v0.20.0-bedpe-results.json。新增入口复用已有导入成功路径测试，未运行testthat。

导入教程formats页在独立R环境执行通过，示例使用随包已观测区间转写，不伪造配对测序记录。
roxygen/Rd、汉语帮助/契约翻译、NEWS、导入教程、API/参数/任务索引及安装版本说明同步。
当前141导出、1573参数、105 Rd主题、990段参考翻译；这些是当前签名与文档清单，不代表全部组合验收。
未执行包级检测、vignette编译、完整图像重绘、整站重建或发布设备矩阵；.site仍是上一版预览，不作为v0.20验收证据。
README未来方向、冻结docs、原始数据和许可证保持不变。

## 首切片时的后续门槛（过程记录）

下一步按版本设计推进真实字符源/位置映射与M0细节图原型审计；公开接口须经原型后确定。
其他VCF/SV、信号、domain/motif适配器按真实来源逐项确定，不一次铺开全部格式。
当前没有字符级比对恢复、语义放大公开入口、Fusion F2/F3或通用字号缩放的实现承诺。

## 本版后续切片完成记录

持续开发授权后完成真实字符与M0的本版原型门槛，运行时提交5a96255：

- read_fasta_lengths追加keep_sequence=FALSE，显式TRUE保留实际去空白字符与大小写；
  默认长度输出、ID、单位声明和文件/物理行来源保持原契约。无字母表/比对/反向互补推断。
- focus在sequence列存在时验证字符串与输入长度，按整数窗口同步裁剪；两次聚焦内容与原始111–130一致。
- 发现M0多层布局投影含缺失文字角度会报NA下标赋值错误；修复纵向和径向同类路径，
  有限角度的可读方向处理不变，缺失角度保留缺失。原失败及修复后的完整M0构建在本地调查记录中核对。
- 隔离字符原型对两个60bp DNA、两个40aa蛋白输出200条真实source-position映射，
  验证DNA反向互补与IUPAC大小写；实测字宽使窄panel选择tiles，不无限缩字。
- M0原型以真实pETDuet-1的101–1900区间重建10 Feature、2引物结合区间、3识别区间，
  原始FASTA和源重叠端点逐条一致；两条装饰leader独立绘制，生物HSP始终为0。
  通用Feature承载线性结合/识别区间，不声称圆图专用Geom支持Collinear。
  已实际查看最终12×8 PNG；10×7只作物理重建检查。三酶各自保留几何，共享名称标签。
- 适配器来源审计和两层字符/对齐方案已形成：ADAPTER-AUDIT-v0.20.0.md、SEMANTIC-RESEARCH-v0.20.0.md。
  PAF确有cg标签，但目前缺同一subject完整字符串；标准蛋白BLAST无gapped列。未虚构差异或域/覆盖度。
- V10-STATUS.md逐项对应长期任务书，保留部分实现、研究与后续门槛，避免把几个版本完成当作全路线完成。

## 最终定向验证与文档范围

check-semantic-detail.R通过真实字符/方向/重复聚焦、物理字宽和M0两尺寸构建；
记录本机原型总耗时（非性能门槛），输出v0.20.0-semantic-results.json、source-letters.csv/PNG、m0-features/annotations/fronts.csv及最终m0-detail.png。
check-bedpe.R在本轮运行时上再次通过。formats独立页面执行通过，包含新keep_sequence示例；
保留可选GInteractions无共同Seqinfo级别的既有warning，不当作导入失败或静默跳过。
991段汉语参考翻译、141导出/1574形参/105 Rd清单、参数图谱/任务索引/双语源页面检查通过。
复用原导入与focus最小测试补必要断言，没有执行testthat。

没有重绘原178张网站图、重建.site/冻结docs、执行浏览器或包级发布流程；
新研究PNG放tools而非冒充正式教程图。运行时跨模块修复的全包/跨设备回归风险留待正式发布验收。
README未来方向、许可证、原始examples/data保持不变。v0.20.0保留专用分支，不自动整合main、push或发布。
