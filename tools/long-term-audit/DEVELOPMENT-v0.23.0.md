# v0.23.0 — M1 多源区间与共线性细节

状态：本版约定范围本地开发完成，未整合、推送或正式发布。
按用户明确批准将v0.22.0 d5dfbbb合入本地main，从main创建v0.23.0。
origin/main仍为69556d8；本次没有推送或触发新的远端CI。

## 任务书与实现

依据v10 §4.10.6 M1：四序列Chord/Collinear主图中选择单区间或真实HSP两端窗口，
进入同一显式Collinear细节面板。magnify_ggchord接受非空多行region，每行独立解析
原始源区间和instance/detail_instance。共享解析与毫米投影不复制为第二套坐标系统。
每个区间各有高亮与两条装饰引线，focus_regions/focus_fronts/leaders中的region_row
对应原始输入行；原始表及同名元数据仍原样保存在focus_region_data。
旧单行调用保持有效，connector='none'不生成引线。

共享Fusion中央面板接受横/纵Collinear，其他受支持面板与既有图层规则不变。
真实HSP由用户显式detail_plot持有，局部裁剪、基因方向和统计归属沿用现有focus流程；
不根据选区伪造比对。没有新依赖、参数改名或旧接口删除。

## 定向验证

check-magnify-pairs.R通过：随包四条真实序列，以ribbon_data_example第一条HSP的两端
作为两个原始源区间，focus后保留一条真实HSP和原有基因链字段。Chord、横向Collinear、
纵向Collinear加原生反射标度的各区间闭边界均正确；两行产生四条独立装饰引线，
不进入跨面板connections/relationship_data。逐行重复实例选择、歧义拒绝、越界拒绝、
非法第二行、单行兼容及connector='none'均直接检查；原始重名元数据保留。
结果v0.23.0-magnify-pairs.json，最终12×8英寸PNG已人工查看。
图中引线可穿过其他序列/几何，密集基因可能重叠；不声明全局避障或所有尺寸质量保证。
已有最小Fusion公开测试改为双区间Collinear成功路径，没有新增成批testthat案例或运行测试套件。

两语Fusion教程含实际成对窗口示例，在独立R环境运行通过；帮助/Rd、汉语参考、NEWS、
支持矩阵、参数清单、安装版本及Rmd源同步。1037参考译段、124页上下文lint、62对源代码和
锚点校验通过；公开签名数量未变。没有编译vignette或重建冻结docs/网站与旧图。

## 剩余范围

本版实现M1的显式独立细节工作流，不等于全部语义放大完成。M2真实字符/CDS phase、
M3多细节面板、自动统计/图层重建、跨原点、Fusion F2/F3、通用Guide/Adaptive仍待后续。
旧远端CI失败记录见CI-69556d8.json：引线回归已在本地main修复，但尚未推送验证；
5个旧lesson engine及178旧图指纹过期待重现。没有完整包检查或发布设备矩阵。
