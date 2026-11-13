# v0.27.0 当前公开入口审查

2026-11-13收尾重查发现教学签名清单没有同步此前已实现的label_segment_layer、
label_repeat、label_segment_colour，audit-public-contracts.R实际退出1；Rd与运行时
本身已有这三项。差异保存在v0.27.0-api-inventory-before.json，没有将旧成功记录
当作当前通过。已只更新受影响入口的清单、汉语缺译、课程重点及自动生成参数表，
再实核运行157绑定/109主题形参通过。当前结果依各R源SHA标识，不等于全部行为
或发布认证；密集布局、当前完整图像及站点重现继续按后续实际记录验收。

这是当前开发签名与职责审查，不是新的1.0冻结或全行为认证。历史v0.16.4冻结候选及
既有发布记录不重写；当前事实分别由实际代码、NEWS、MIGRATION及版本记录约束。

## 实际核对

`audit-public-contracts.R` 加载当前源码，逐一比较157个导出绑定与教学清单。其中135个
由ggchord拥有，其余22个保持上游原绑定，包括PositionIdentity对象。函数默认表达式
按R语法树解析再比较，不运行表达式，不把deparse换行空格当破坏性改变。当前规范化
签名差异为0。代码来源由结果文件逐源SHA256确定，HEAD只提供检查时提交背景。

对函数对象的相等检查只用于确认上游绑定未换实现，不能据此断言语义相同；例如
after_stat与after_scale具有不同staged语义。源声明扫描未发现把参考软件/机构品牌
写入函数或参数名；该扫描不是一般文本、全部动态调用或许可认证。

## 当前保留决策与职责

| 入口组 | 职责与保留决定 | 实际依据/验证边界 |
| --- | --- | --- |
| ggchord及逐层data/mapping | 构造默认输入；每层保留独立来源和身份。第三位置仍为gene_data，与feature_data互斥 | 实际签名、annotation-input.R、MIGRATION统一注释；签名审计不等于任意多层组合认证 |
| gene/feature及plasmid快捷入口 | gene为有方向语义的快捷层，feature容纳未知链和多种形状；plasmid保留不同默认样式 | api-contract.R及已有真实目录检查，不为减函数数合并改变默认表现 |
| geom_track与bar/col/line/tile | 通用层拥有受控Geom/Stat桥接；快捷入口保持对应默认。bar为count、col为identity | TRACK-GRAMMAR记录、track-bar-migration.csv及估计统计定向证据；不承诺任意原生透传 |
| 六种coord | 区分投影、可用实体与序列语义，不以统一名字隐藏限制 | 公开coord帮助与既有视图证据；不把Circular的禁止关系层视为待自动放宽 |
| scale/guide与颜色双拼写 | 规则与guide由scale拥有，用户scale优先；colour/color拼写保留 | 公开帮助及主题/guide工作包；保留ribbon与link各自几何语义，不因近似名称删除 |
| theme及五个预设 | 保留109形参，便捷参数与注册元素走同一消费点；后续theme逐字段覆盖 | THEME-DESIGN和check-theme-system记录；精简签名暂无获确认的新删除计划 |
| 适配与来源工具 | reader负责来源格式和明确坐标转换，转换/清理/筛选/聚合分开 | 各适配工作包真实来源记录；不以reader数量或临时可选库存在宣称任意文件支持 |
| focus/segment/project/position/映射 | 区间、实例与坐标映射各有目的，不按前缀或内容相同猜来源 | 已有复审修复、真实映射证据；未知长度/方向/支持度不填造值 |
| Fusion/magnify/字符注册/CDS | 组合容器、显式关系、焦点与字符映射分开，保留嵌套身份和原记录 | F2/F3/M2/M3定向记录；复杂高密度及全面物理质量仍待补证 |
| get_chord_layout/export/view/diagnose | 完整布局、稳定导出表、真实尺寸预览与诊断目的不同，双布局getter不合并 | layout-contract及output教程；诊断不修改来源或保证全部约束，正常保存仍ggsave |

## 兼容边界

- `gene_data` 仍支持，包括第三位置；仅构造时发迁移消息，无已确认删除版本。
  新例优先feature_data；typed gene/CDS选择与显式gene层见MIGRATION，不同时传两个构造入口。
- 已删除的旧形状、offset及标签参数继续报替代方向，不恢复静默兼容。既有删除版本及
  原minus链offset迁移取负写法以MIGRATION和api-contract为准。
- theme快捷形参、常用geom快捷层、拼写别名及get_chord_layout继续保留。没有仅为
  缩短参数表新增一组碎入口，也没有在本次审计中宣布新的弃用/删除期限。
- optional依赖保留明确错误/跳过语义；签名和文档相符不证明所有依赖版本行为一致。

## 未关闭

全部公开行为的跨模块、不同依赖与设备认证尚未完成。后续aaa730c运行时的50课/
178单PNG已经真正重现且源检查通过，新R最终HTML/归档/链接及4个代表浏览器状态亦已通过所列范围；这不重新认证本篇
签名审查以外的全部语义。制图质量、复杂Fusion/Adaptive仍有未验证项。
进入1.0的最终纳入范围及新删除计划没有在此擅自决定。包级及完整设备/平台检测仅在
用户明确授权GitHub正式发布后执行，本次审查没有运行这些流程。
