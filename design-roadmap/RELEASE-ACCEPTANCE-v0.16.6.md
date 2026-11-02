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
