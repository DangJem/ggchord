# 首页上下布局与 Radial/Sankey 示例优化

## 授权与范围

基线为 v0.27.0/3a6b412。用户要求八项首页视觉调整及分析优化 Radial/Sankey
旗舰。后续答复明确继续仅做源码、构建与图片检查，不操作浏览器；桌面、窄屏、
短窗口和放大状态的响应式源规则可以修改，但实际视口/缩放/computed样式未验收。
继续原 checkout，本地阶段提交，不整合 main、push、部署或正式发布，不跑包级
检查、vignette重建或设备矩阵；原始数据、R包实现与冻结docs不改。

## 问题分析

1. 旧首页左品牌栏占约25%宽度，右四模块再分两列；多序列模块内部四张图再分
   2×2，单图实际宽高都很小。原 `23svh` 图片区域进一步随短窗口压缩。
2. 品牌区与模块整体垂直居中，介绍/输入要求/预览争夺首屏；图像内容与用途主次
   不清。后期覆盖虽移除部分固定height，原svh最小高度和重复媒体规则仍存在。
3. 两张派生旗舰沿用原Chord的100mm高左侧Identity标尺、底部重复Sequence图例及
   12×9画布。Sankey四条序列各在左右列纵向堆叠，ratio保持比例时主体窄；Radial
   的四正交辐条和旋转ID进一步占边缘空间。字体/箭头缩到网页后难以区分。
4. 同一完整科研图不能直接用统一裁切解决缩略图问题：会截断来源区间/文字。
   本轮分别生成完整图和首页预览，并明确预览省略的文本/图例，不改科学数据。

## 实现

- `tools/site_components.py` 生成真实首页结构；`website/home.json` 保存英汉同结构
  文案，`home-previews.json` 保存预览与对应完整图。没有单独造一套替代发布页面。
- 介绍区左文右Logo，Logo224px（中宽度200px）；用途作唯一H1，删除Logo下面的重复
  名称。内容最大1280px，上方介绍/三项核心能力/下方全宽模块，不设视口高度。
- 模块名称→用途→预览→图型名→进入教程。桌面两列，800px以下单列；560px以下
  四视图也逐图单列，保证小宽度下不会继续把四图挤成极小缩略图。所有图片保持
  导出比例、height:auto与contain，无裁切、固定svh或整体居中。输入字段留模块页。
- 保留Logo、四模块、“灵活的参数格式”、暖白/青绿与科研图原色。模块三的面向
  用户入口改为“序列比对点图”，保留Dotplot；技术HSP字段/统计语义不重命名。
- 首页按钮“进入教程”为用户本轮指定文案，其他说明书口吻保持。预览图均可点击
  打开对应完整图；同模块有实际用法入口和共享预览源码，不将别的图当同一完整例。
- 完整Radial：11×10in/144dpi，45°起始、水平序列名、较厚CDS，20条原CDS与12条
  原HSP保留。名称输入按每序列优先非假设蛋白、再按长度选一条；布局策略仍可能
  隐藏拥挤文字，不保证每个候选名都显示。
- 完整Sankey：10×12in/144dpi、显式列间距5.2，16条原CDS与10条已选HSP保留。
  名称仍每序列一个代表，几何和名称字号增大。左右源副本不增加观测。
- 两图Identity移到下方130×4mm横向标尺，保持0/50/80/90/95/100及数值区间比例。
  序列ID已直接标出，移除重复Sequence图例；Strand保留。没更改包的默认外观。
- Chord/Collinear完整旗舰及Dotplot/Fusion原完整图、13完整参考图继续原字节/来源。
  修改site-flagships.R仅涉及flagshipmulti片段；保留图逐片段等价，未改旧生成hash
  冒称新图。`check-manual-assets.py`分别核对新渲染和旧来源。

## 首页预览数据与显示选择

| 预览 | 数据来源及保留 | 仅显示层的调整 | 完整图 |
| --- | --- | --- | --- |
| Chord | 原Chord来源12 HSP/20 CDS及原序列曲率/半径 | 不绘CDS名称/图例/密集坐标字；适配线宽 | 原Chord旗舰 |
| Collinear | 原10 HSP/16 CDS、原跨行路由与显示方向 | 不绘CDS名称/图例/密集坐标字 | 原Collinear旗舰 |
| Radial | 原12 HSP/20 CDS及显示方向 | 不绘CDS名称/图例/密集坐标字；较厚箭头 | 新Radial完整图 |
| Sankey | 原10 HSP/16 CDS及左右副本 | 不绘CDS名称/图例/密集坐标字 | 新Sankey完整图 |
| Circular | 现有真实pUC19全部9条识别特征 | 名称仅CDS/rep_origin，保留源feature_color/segments；不添加实验信息 | pUC19完整参考图 |
| Dotplot | 原四序列全长及全部32 HSP | 轴段中心仅标原ID；完整坐标在完整图；线宽调整 | 原32 HSP矩阵 |
| Fusion | 现有双面板示例的全长表及first fusion_links真实HSP | 中心名称和Query/Subject显示名、预览装饰 | 原coursefusion双面板图 |

不补对称HSP、不改链/坐标、不裁切PNG。Fusion放弃四面板压缩候选，其径向面板
被挤成小点且曾未收敛；当前用现有双面板来源，最终渲染无该警告。候选过程PNG
只在临时目录，不纳入仓库。

## 验证与失败事实

- `Rscript --vanilla tools/site-audit/render-home.R`：实际生成2完整+7首页PNG；通过
  公开布局的layer_inputs逐项比较前后HSP端点/pident与CDS accver/start/end/strand/anno。
- 首轮在homecollinear遇到未分型feature_data错误：改为给geom_gene显式原CDS表，
  不伪造类型或改原表。原partial记录为home-rendering-initial.json，非成功记录。
- 实际图片已检查：完整Radial/Sankey、7首页预览；修复预览Sankey坐标文字拥挤、
  Dotplot轴标签边缘过长、Fusion竖列名称相碰。不称图像查看等同浏览器查看。
- `Rscript --vanilla tools/check-documentation-lessons.R --page=multi-sequence-views
  --update-outputs`：更新页面片段构造与两种导出尺寸说明，准备通过；不跑包测试。
- `python3 tools/check-documentation.py --source-only`：93语言对、157导出/121帮助、
  代码配对与既有资源核对通过。
- `python3 tools/check-site-design.py --write`：76受控色彩组合与6静态条件通过。
- `python3 tools/site-audit/check-manual-assets.py`：20完整资产、7首页预览、原始来源
  与显式尺寸核对通过。
- 初次CDS来源比较取错component键，gene空表被误认为相同。已拒绝该证据并改用
  gene_poly、强制非空；`render-home.R --verify-sources-only`实际核对Radial12HSP/20CDS、
  Sankey10HSP/16CDS全部源字段相同。原PNG生成hash保持，旧空检查与补正分别记录。
- 页面表格输出刷新带来一处无关null→空对象序列化，已撤回至基线。进行中的HTML
  构建因此被快照guard正确拒绝，退出1；记录home-build-initial.json，不计通过。
- 当前本地原生双语HTML、生成链接/下载与旧锚点核查待完成，不预写通过。

## 明确未验收

用户选择不操作浏览器。本轮未检查实际桌面/窄窗口/短窗口/200%缩放、键盘/主题
交互与computed样式。CSS断点和自然高度的源码检查不代表这些状态已通过。
没有正式发布、部署或包级设备/平台/字体认证；本轮不改底层几何算法或原始记录。
