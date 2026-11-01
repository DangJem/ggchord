# 测试组织与WORDLIST调研

2026-11-02；当前开发分支 `codex/v0.16.6`。用户要求调研，不将此解释为删除测试或
新增图形快照的授权。没有执行本包或外部项目的测试、build/check和设备矩阵。
本地main已恢复至v0.16.5的 `d697345`；今后整合须用户明确批准。

## 外部项目观察

从各项目官方GitHub默认分支先固定commit，再读取完整、未截断的递归tree与代表
文件；不是CRAN发布包对比。机器可核对的commit和计数见
`audit-generated/test-survey-v0.16.6.json`。test文件同时匹配test-与test_；数量表示
源文件，不能解释为case数、执行通过数、覆盖率或质量排名。这里只概述，不复制代码。

| 项目 | test文件 / SVG快照 | 实际组织与可借鉴点 |
| --- | --- | --- |
| [ggplot2](https://github.com/tidyverse/ggplot2/tree/b638d654c8189db82a62307fefa84b1da1819da5/tests) | 142 / 314，另99份文本快照 | 按aes/geom/coord/scale等模块分文件，数值行为和错误文本有直接断言；helper-vdiffr独立处理可选视觉依赖。适合借鉴分工，不宜照抄其规模 |
| [gggenomes](https://github.com/thackl/gggenomes/tree/752e0e19da04b6813f0524d784843901ac64dcb2/tests) | 2 / 0 | test-read_seq_len覆盖FASTA/GFF/GenBank长度读取，test_vars覆盖变量检查；真实小输入集中在tests/data。当前公开测试范围窄，不能因其少就推导ggchord只需两份或无需几何回归 |
| [gggenes](https://github.com/wilkox/gggenes/tree/125b47a1a9cb0071c2053c810ddc01a80056fa42/tests) | 8 / 16 | 按gene_arrow、subgene_arrow、feature、terminator等职责分文件；helpers复用基础图、空设备绘制及几何检查，兼有smoke和vdiffr。另有独立spelling.R，启用可选拼写检查 |
| [ggrepel](https://github.com/slowkow/ggrepel/tree/4a3100bf3970b0aadb549326bb1925c90d9d8dac/tests) | 8 / 16 | 按seed、side、theme、grob顺序等行为组织；setup重置主题；兼有图形快照、绘制及连接几何数值断言。文件少也可以包含多设备调用和内部细节，不自动等于轻量 |
| [patchwork](https://github.com/thomasp85/patchwork/tree/6b1d88ce1da1c5cae3818d984edf80dc0bb3de8c/tests) | 3 / 71 | 按arithmetic、collect、layout三个大主题，复用基础图和vdiffr helper；一个layout文件覆盖很多可视组合。少文件与少测试不是同一件事 |
| [ggraph](https://github.com/thomasp85/ggraph/tree/6c4ce81e460c50a16f9cd97e0b3a089f36901316) | 当前tree无tests目录 | 存在R-CMD-check/coverage工作流，但该快照没有tests；不能由此证明质量充分或完全不做验证，也不建议仿照删除目录 |
| [ggforce](https://github.com/thomasp85/ggforce/tree/128001400bee8ae1a27c772b6d9677660db01778) | 当前tree无tests目录 | 同样区分公开tree事实与质量判断；其不存在tests不是ggchord删除严重缺陷回归的依据 |

代表源码：[ggplot2视觉helper](https://github.com/tidyverse/ggplot2/blob/b638d654c8189db82a62307fefa84b1da1819da5/tests/testthat/helper-vdiffr.R)、
[gggenes基础图/绘制helper](https://github.com/wilkox/gggenes/blob/125b47a1a9cb0071c2053c810ddc01a80056fa42/tests/testthat/helpers.R)、
[ggrepel侧别测试](https://github.com/slowkow/ggrepel/blob/4a3100bf3970b0aadb549326bb1925c90d9d8dac/tests/testthat/test-side.R)、
[patchwork布局测试](https://github.com/thomasp85/patchwork/blob/6b1d88ce1da1c5cae3818d984edf80dc0bb3de8c/tests/testthat/test-layout.R)。

testthat官方建议测试文件与源码对应，并明确helper/setup/_snaps的不同用途、避免
测试文件相互依赖。这是组织建议；ggchord有大量私有布局文件，不需要机械新增同等
数量的test文件。[官方特殊文件文档](https://testthat.r-lib.org/articles/special-files.html)

## ggchord当前问题

静态读取：10个test文件、1个13行helper、1个入口；testthat目录共719行，45个
test_that代码块。tests总文件字节52,670，无快照/二进制图片。没有执行这些case；
循环中的运行次数及expect数量不能由静态块数推断。总体不算大，但分类存在交叉：

- plot-api同时含三标签策略、region/feature/link及来源映射；feature-engine仅21行，
  link-api仅28行，track-api仅14行，职责分散。
- radial-geometry同时含设备所有权、Unicode文字度量、标签几何及圆图混合注释；
  其中设备/cache又与layout-ownership交叉。文件名不能完整反映内容。
- data-api为230行、15个case；第一个case同时检查validate、clean报告与constructor
  错误，preparation/loci案例串起多项工作流。失败时需要读较长上下文才能定位合同。
- 基础成功路径、严重回归与有限几何检测混在一起，缺少简短的来源/缺陷说明；
  F52/F53/F54/F55等应能定位，但不能单凭“老测试”或内部helper调用就删除。
- 默认helper只做小数据与一次build，符合轻量方向。现有内部碰撞检测可验证公开
  布局的方向/保留结果；精确浮点坐标、grob结构或所有等价参数组合不宜作为新门槛。

## 建议（待确认，尚未重组或删改tests）

优先按五个公开责任域收拢，保留标准testthat入口和小helper：

```text
tests/testthat.R
tests/testthat/helper-fixtures.R
tests/testthat/test-data.R         # validate/clean/import/utilities/workflow/projection
tests/testthat/test-grammar.R      # constructor/core geoms/mapping/scales/axes
tests/testthat/test-annotations.R  # feature/labels/primer/restriction/direction
tests/testthat/test-collinear.R    # anchors/segments/local reversal
tests/testthat/test-output.R       # export/cache/device ownership/physical track evidence
```

这是一种可审阅的组织草案，五不是质量配额。实施时第一步只移动和归类、保留45个
已有case及断言，再单独判断重复。每个已确认严重缺陷保留一个最小公开反例，来源
可用短注释关联F编号；当前不能只凭相似标题合并“复用layer”和“派生plot”两种所有权。
数据格式各自的一个成功路径也不能相互替代。未知/空ID的多种政策可审查代表性，
但在来源和严重性尚未核对前不删掉边界。

继续维持AGENTS的既定政策：tests只保留最小公开行为与严重回归；性能、复杂布局、
完整设备验收仅在tools显式运行，不机械搬运被删测试。不要为减少文件数而把所有
测试塞进一个长文件，也不因外部项目使用vdiffr就新增一套快照维护。
外部包针对每参数/内部grob/精确数值的做法与本项目现行限制不同，不直接照搬。
后续若批准重组，普通开发只做parse及移动前后case/断言清单核对，完整testthat仍
须符合既定发布授权条件。

## WORDLIST结论与维护

它有用途：这是文档拼写检查的自定义词典，不影响图形运行、数据导入或标签布局。
官方规定UTF-8、每行一个词；工具会把这些词加入字典，自动update_wordlist会将
候选加入词表，因此必须人工区分正确术语与真正错误。
[WORDLIST文档](https://docs.ropensci.org/spelling/reference/wordlist.html)、
[包拼写检查文档](https://docs.ropensci.org/spelling/reference/spell_check_package.html)。

原文件只有BLAST/ggchord/ggplot三项，并混用CRLF且末行无换行。仓库没有显式
spelling测试/CI入口，也未声明spelling依赖；保留这个文件不会自行触发检查。

本轮将其更新为141个逐项审阅的当前术语/缩写/正确专名及文档已有的合法英式拼写，
统一UTF-8/LF，保留旧3项。词来自当前英语Rd、DESCRIPTION、README、迁移和教程，
例如accver/qaccver/saccver、CDS/GFF/PAF、Biopython/Biostrings、ggproto、collinear、
colour/centre。没有把纯汉语词或未审阅候选批量加入，不改作者、原始记录或显示标签。

新增 `Rscript tools/check-spelling.R` 显式入口：仅阅读Rd/元数据及英语文档，使用
现有可选spelling包，不安装依赖、不写文件、不运行tests/knit/build、不接入CI或
增加Suggests。英语checker排除纯汉语导航文字，不能据此声称检查了汉语拼写。
DESCRIPTION目前无Language字段，现有工具默认en-US；本轮没有借机改变元数据语种。

针对性运行后仍报告7个Rd候选和2个教程候选（canonicalized、oversized、overwide、
pre、spacings、undisplayed、unexpanded、projectable、unshown）。逐项对照为正确
上下文/合成词或pre-computed的拆词；没有确认真实拼写错误，不以填满白名单强求0。
输出位于忽略目录 `.codex-tmp/test-research-v0.16.6/spelling-review.txt`。该结论仅
覆盖本轮选定文档，不冒充全仓库或网站的拼写验收。
