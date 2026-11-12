# GFA/ODGI 探索性桥接结论

对应v10 §6.6：探索节点/路径到已有视图的子集转换，不是完整组装图Coord承诺。
本阶段没有新增公开reader，未调用比对器/组装器/ODGI或生成新的科学关系。

## 原型、输入和复现

实际使用作者固定v0.8.4提交3f59d596cff3bbe9ddff3862e3212812254496e5的原始
人类LPA图（13条contig路径、7个单倍体组装的局部区域；具体来源以原path为准）。
作者教程说明来源：<https://pangenome.github.io/odgi.github.io/rst/tutorials/presence_absence_variants.html>。
原始URL、字节SHA256及MIT归属见graph-path-results，原许可完整保存在
GRAPH-PATH-LICENSE.txt。没有复制上游实现或重新分发完整原图；CSV/PNG只是明确
选择的源路径步骤展示，保留原始ID/顺序/出现索引，不改源字符或节点。

明确的开发复现命令（Python环境需自行具备可选gfapy1.2.3；本次用/tmp隔离环境）：

```sh
python3 data-raw/prepare-graph-prototype.py /tmp/ggchord-graph-cache
/path/to/optional/environment/bin/python tools/long-term-audit/prototype-graph-paths.py /tmp/ggchord-graph-cache/LPA.gfa
Rscript tools/long-term-audit/check-graph-paths.R
```

实际本轮使用重新下载、匹配哈希的/tmp/ggchord-odgi-LPA.gfa；gfapy只通过公开
Gfa.from_file接口解码。3751节点/5195连接均为0M、13路径共202806个源步骤。
原型逐转移复核真实有向边或其反向互补遍历，逐步骤核对真实节点字符长度；
累计得到的是原始路径全长，不是从一个局部端点推断整条染色体长度。路径ID内
的双下划线及contig名称保持opaque，不猜其参考组装坐标或自动改名。

两个明确命名路径只展示各前100个原步骤，200行CSV中保留node_id/source_step/
node_occurrence/path strand与累计一基闭合区间；全13路径长度另表保存。真实
seq/feature/focus公开API构建通过，PNG已实际查看。两条所选前缀实际均为+遍历，
不为展示正负颜色伪造反向步骤；所有原路径的真实方向在原型中逐条消费。

## 为什么没有直接发布read_gfa

- 当前数据的0M路径可复用现有Seq/Feature；同节点多次出现，node ID不是feature
  实例。原数据每路径存在773–887个重复节点ID，因此自动按节点做pair连接会
  产生大量歧义和错误科学暗示。原型生成0条科学pair链接。
- GFA规范允许非零/未知CIGAR重叠、无字符S、GFA1.1 W、1.2 J和GFA2；不能简单
  累加端点得到可靠路径坐标或推导identity/支持度。本原型遇非0M/缺字符明确失败，
  没有宣称是通用GFA解码器。当前gfapy版本对W/J范围也未在这里验证。
- 原生ODGI二进制需成熟工具先导出有来源的路径/GFA子集；本轮没有安装或认证
  其二进制，不把导出命令文档冒充实测。官方view/paths命令是可选上游桥接路线。
- 决策：探索门槛已有真实原型及限制，优先使用显式路径/Feature子集和已验证的
  普通表接口，不增加空壳reader或Graph Coord。未来如公开通用入口，须补重叠/
  缺长/W/J/图身份/显式occurrence配对语义及真实多格式输入，不自动成为版本承诺。

规范依据：<https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md>；
成熟解码接口：<https://gfapy.readthedocs.io/en/latest/tutorial/gfa.html>；
ODGI视图转换：<https://pangenome.github.io/odgi.github.io/rst/commands/odgi_view.html>。
这些是格式/公开接口借鉴，函数/参数保持通用命名，不复制第三方代码。
