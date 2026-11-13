## 真实比例的径向辐条 {#radial}

随包四序列数据保留原始长度和CDS链。以下显式展示筛选选择长度至少900 bp、identity
至少75%的已观测HSP，并取每个有向序列对最长两条。计数轨道统计所选query起点观测数，
不将其解释为深度或覆盖度。

{{code:radialdata}}

{{code:radialview}}

{{figure:viewradial}}

辐条顺序来自`geom_seq(seq_order=...)`。默认按真实长度比例展示；只有显式
`length_mode="equal"`才归一显示长度，源坐标保持不变。`start_angle`从正x轴逆时针计角，
`clockwise`控制辐条顺序，不改生物链。特征和轨道的正偏移指向辐条逆时针法向。
物理Track宽度和源统计复用共享band与Stat系统。

首尾序列只有在**存在HSP**时才在闭合通道连接；没有该关系时通道保持空白。
N=2每条原始HSP只绘一次，N=1报错。反向HSP保留完整、有顺序的前缘及扭转。
新Radial采用闭区间边界`[min-1,max]`，所以单位置HSP仍有宽度；轴位置仍是一基。
已有视图继续采用各自文档约定。中心孔只是显示留白，不删除序列；序列名位于辐条外端。

{{prepare:coursecomparison}}

{{code:courseradialangle}}

![仅改变首辐条角度](../assets/courseradialangle.png)


Ribbon可用显式`pairs`或`adjacent_only=FALSE`选择非邻接关系；非邻接点连接须指定
`between`。命名A–B–A实例沿用[明确路由契约](foundation.html#instances)。
用户控制点使用最终x/y布局单位。端点渐变和subject尖头仍须显式选择，默认保留Identity填充。
## 导出与反向映射 {#export}

{{code:radialmap}}

导出保留`radial_frames`、原始源区间与`radial_endpoints`中的完整投影前缘，以及连接
筛选、实例路由、轨道统计和标签可见性。Radial反查必须给出所选骨架上的x、y，避免
仅凭多条辐条共享的单个x值猜位置。focus/segment输入保留源原点；径向显示不将多个
窗口拼为一条辐条。

## 练习与支持边界 {#support}

只改 start_angle 不改其他参数，核对源坐标及方向。N=1、shared branches、单圆专属层和多窗口拼接为一条辐条不支持；逆映射要给骨架上的 x/y。
