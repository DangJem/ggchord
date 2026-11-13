搜索切点与选择显示集合是两步。用相同三个 pUC19 命中比较引线、文字内容、顺序和内外侧，不修改 DNA 或切点。

## 保留真实切点与搜索条件 {#input}

本页独立运行。显式使用 EcoRI、BamHI、HindIII 和 methylation="none"；这是搜索条件，不是样品甲基化测量。enzyme_site_count 针对完整输入序列，筛选显示行后仍保留该统计。

{{code:comparefinishdnabase}}
{{code:comparefinishsites}}
{{code:comparefinishsiterows}}
{{output:comparefinishsiterows}}


## 同一切点比较引线 {#restriction-labels}

radial 默认根据文字框与圆周顺序排布，必要时使用短径向段与扇形引线；elbow 用一个折点，straight 连接切点短线端与最终文字边缘。三图只改 leader，位置和标签分组相同。

{{code:comparerestrictionradial}}
{{figure:comparerestrictionradial}}
{{code:comparerestrictionelbow}}
{{figure:comparerestrictionelbow}}
{{code:comparerestrictionstraight}}
{{figure:comparerestrictionstraight}}


## 内容、顺序与侧向分别改变 {#text-content}

以下每幅都回到 radial 基线，分别只改 label_style、label_order 或 label_side。酶名加粗默认表示它在完整输入序列中只切一次；显式 fontface 会替代该自动区别。位置来自实际切割结果，不是识别区间中点。

{{code:comparerestrictionenzyme}}
{{figure:comparerestrictionenzyme}}
{{code:comparerestrictionposition}}
{{figure:comparerestrictionposition}}
{{code:comparerestrictionorder}}
{{figure:comparerestrictionorder}}
{{code:comparerestrictioninside}}
{{figure:comparerestrictioninside}}


## 关键参数与筛选含义 {#parameters}

| 参数 | 默认 | 作用 |
| --- | --- | --- |
| label_style | enzyme_position | enzyme_position/enzyme/position |
| label_order | auto | auto/enzyme_position/position_enzyme；自动随侧向镜像 |
| label_side | outside | outside/inside/auto |
| leader | radial | radial/elbow/straight；不改变切点 |
| min_label_gap | NULL | 按设备文字框计算；显式值为基因组比例间隔 |
| tick_length | 0.04 | 局部法线方向短线长度 |
| label_offset | 0.18 | 文字局部法线位移 |
| linewidth / label_size | 0.28 / 3.2 | 线宽和文字字号，含义不同 |
| fontface / family | NULL / NULL | 固定字体设置；显式 fontface 取消自动唯一切点加粗 |

find_restriction_sites 返回全部匹配和切点，filter_restriction_sites 选择显示集合。unique、six_plus、unique_6plus 与 reference 有各自条件；reference 仅在精确匹配的参考上使用保存集合。methylation=auto 只对已审计的精确内置参考使用背景，未知样品需明确假设。



## 练习、密集图与限制 {#adapt}

只增大 min_label_gap，观察标签与引线是否向外扩展；再回到默认，用 label=FALSE 只显示切点。保持完整搜索结果，另存显示子集。小图放不下时增加实际导出尺寸，而不是删除源记录。组合特征与引物后重新检查拥挤，见[密集 pETDuet-1](circular-dense-case.html#complex)。计算切点不是实际酶切实验的验证。
