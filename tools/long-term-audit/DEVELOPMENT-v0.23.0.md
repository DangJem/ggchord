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

## 2026-11-12 审计修复

用户要求先修正v10复审结论，再继续开发。本版维持0.23.0.9000，修复以下已支持行为：

- 显式Track group优先于离散样式推断，不再把同一统计组拆开。
- Bin的来源成员在原生Stat阶段确定并随结果保留，包含外边界和数值容差；
  pad/drop和数据Position不再重新推断成员。用户breaks/binwidth回调每组仍只求值一次。
- 普通Chord/Collinear和已有Radial/Sankey共享过滤后检查；只有实际显示的声明HSP参与
  重复诊断。bundle结果通过已有成员报告回查原始观测；重复显示仍要求各参与层显式同意。
  instance_routing增加retained，repeated_display只描述实际显示的重复来源。
- magnify的可选instance/detail_instance字段精确访问，普通同前缀元数据不再选中副本。

审计原始结果仍保留；修复复现命令：

```sh
Rscript tools/long-term-audit/review-v10-2026-11-12.R tools/long-term-audit/review-v10-2026-11-12-fix-results.json --check
Rscript tools/long-term-audit/check-magnify-pairs.R
```

四项defect_observed均为FALSE。额外定向核对了pair/邻接过滤、显式允许重复、bundle来源，
以及左右闭分箱的边界容差、pad/drop、数据Position和单一grob构建。
M1原有脚本改为逐源行核对strand/原始端点，增加真实HSP裁剪后两端的同一进度与origin检查。
物理锚点独立使用实际骨架插值和实测panel边界核对，覆盖Chord、横/纵Collinear及原生反射；
没有冻结精确设备浮点坐标。单张PNG重现并人工查看，图像内容未变；不声称全局引线避障。

复用既有最小Track与Fusion测试，补入分箱来源和声明HSP重复的单一回归；未运行testthat。
NEWS、roxygen/Rd与当前状态同步，git diff --check通过。未运行包级检查、发布矩阵或整站构建。

分箱适配器局部调用ggplot2未导出的compute_bins/bin_cut以共用原生归属，不复制分箱算法。
已核对[ggplot2 4.0.0源码](https://github.com/tidyverse/ggplot2/tree/v4.0.0/R)，
本机4.0.3定向验证通过；将来升级ggplot2时须复核这一适配点，不能宣称任意未来版本兼容。
无声明ID的跨层来源继承仍属后续设计，未将等值数据擅自去重。

## 后续开发交接

按v10 §6.3.1与A.7，下一阶段优先补齐真实画布下的Guide联合布局，再推进可选Adaptive。
先复用guides.R的测量、render context与既有gtable重测机制，设计Guide约束和候选原型，
覆盖同侧多组图例、长标签、Identity比例色条，以及用户显式行列/方向/位置优先。
达到可行性门槛后再确定公开入口；fixed基线、字体下限、失败诊断和Fusion子面板测量是边界。
M2/M3与高级Fusion继续保留，不将它们与本轮4项修复混为完成。下一版实现须先获准整合本版，
再从main建立对应版本分支和完整设计文档；本次尚未开始下一版运行时开发。
