# v0.20.0 Semantic Enrichment 开发记录

状态：开发中，首个BEDPE适配切片已实现并完成定向验证，不表示整个Semantic Enrichment阶段完成。
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

## 后续门槛

下一步按版本设计推进真实字符源/位置映射与M0细节图原型审计；公开接口须经原型后确定。
其他VCF/SV、信号、domain/motif适配器按真实来源逐项确定，不一次铺开全部格式。
当前没有字符级比对恢复、语义放大公开入口、Fusion F2/F3或通用字号缩放的实现承诺。
