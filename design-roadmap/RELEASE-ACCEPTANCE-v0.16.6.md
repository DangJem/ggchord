# v0.16.6 正式发布验收

2026-11-02，用户明确授权整合main与GitHub正式发布；当前发布候选a3ffe41。

- 完整devtools::test：166 PASS，0 FAIL/WARN/SKIP。
- 标准R CMD build与R CMD check：Status OK；完整Suggests安装后检查，无跳过依赖的环境覆盖。两份vignette重建、帮助例子、包内测试及PDF手册通过。
- 9个来源案例×PNG/PDF/SVG，共27文件通过，有SHA-256、单页PDF和SVG XML记录。覆盖三坐标、auto/radial/arc、定量轨道、两个独立gene/ribbon层和变化几何。准备输入的倒向HSP诊断不是设备警告。
- Unicode圆图采用Quartz PDF，其余标准PDF；本机不可用的Cairo不宣称通过。
- 39对网站页面、47个当前PNG、13个独立质粒案例；865译文、130候选导出与93个Rd主题。文档源检查通过，完整产物及远端结果继续补录。
- 声明下限R>=4.1、ggplot2>=4.0不等于本次实际测试组合；密集文字可隐藏/跨图层交叠，无全局零碰撞承诺。组合视图与测试精简仍为待确认构想。
- 版本分支按用户新要求去掉codex/前缀；15个本地、7个既有远端版本分支已改名，当前v0.16.6已推送。历史用途后缀、备份分支、发布标签不改。

远端检查、main整合、Pages部署与正式Release尚待完成；不以本地通过冒充已发布。

## 远端平台结果

[四平台检查](https://github.com/DangJem/ggchord/actions/runs/41414761260)全部成功：
macOS release、Windows release、Ubuntu release与Ubuntu devel。对应R包候选a3ffe41。
29个文档准备环境通过（已安装Bioconductor可选依赖）；原始倒向HSP诊断保留。
源码包SHA-256：`972b9aec161bf9c527e1569e68e46ac20a2d92383151f5ee2261c6b3f1860156`。

## 正式Release

[Release](https://github.com/DangJem/ggchord/releases/tag/v0.16.6)于2026-11-02
16:48:23（Asia/Shanghai）正式发布。标签指向58f9ed9，双语说明现代汉语在前、English在后，
非草稿、非预发布；源码包及SHA256SUMS.txt服务端摘要与本地文件一致。main已快进整合推送。
[配套网站部署](https://github.com/DangJem/ggchord/actions/runs/41415616547)正在执行。

正式Release两附件已下载回读，源码包SHA-256与SHA256SUMS.txt一致。

## 配套网站最终结果

[预发布网站验收](https://github.com/DangJem/ggchord/actions/runs/41415101818)成功。
[正式Pages部署](https://github.com/DangJem/ggchord/actions/runs/41415616547)的source-check、
preview（完整构建/产物检查）、deploy全部成功，来源58f9ed9。
线上英语/汉语首页回读确认0.16.6、三张展示图及发布版标识。
本地39对页面、865译文、130导出/93Rd与链接/锚点/搜索完整检查通过。
后续仅发布记录与历史索引提交，不改变已验收R包、vignette或网站输入，发布标签不移动。

## CRAN修正后的GitHub同版本重发（2026-11-02）

用户明确授权替换当前Release。修复提交1741fcf及授权记录d49d6aa已整合main并推送，
v0.16.6标签从58f9ed9更新为d49d6aaf188d2658c7f2537e9bc7d06011a15a61。
[原Release页面](https://github.com/DangJem/ggchord/releases/tag/v0.16.6)保留，
双语重发说明均已补充修复、验证范围、安装、重新下载及未获CRAN接受的边界。
源码包和SHA256SUMS.txt均已替换，当前包4,834,097字节。
SHA-256：`0a5aaa55b589e86fa865694a8d065fd66a2a8e638b47e93e9e021f968e7195a7`。
远端tag、完整双语说明、非草稿/非预发布状态、服务端摘要已回读；两附件已重新下载，
源码包SHA-256与下载的校验文件一致。旧版事实及旧摘要保留在前文，不覆盖历史。
运行时代码和公开签名不变，47张实际重建PNG逐字节一致，复用首次27设备结果；
修正候选本地as-cran、双语vignette、完整依赖/测试及最终资源闭包证据见
[CRAN修正记录](CRAN-FIXES-v0.16.6.md)。最后仅审计记录提交不移动新标签。
本次未向CRAN webform提交，也未重新部署网站；GitHub推送触发的CI结果不冒充已有完成结果。
