# v0.16.5 正式发布验收

2026-11-02；用户已明确授权整合 main、GitHub 正式发布及新网站同时上线。

- `devtools::test()`：162 PASS，0 FAIL/WARN/SKIP，包括稀疏引物/酶切混合引线实际相交回归。
- 标准 `R CMD build` 与 `R CMD check`：Status OK，完整安装帮助例子、包内测试、双语 vignette 重建及 PDF 手册成功。首轮发现的帮助示例错误、缺失参数说明与未限定 utils 调用均已修正；未跳过 Suggests。
- 9 个有来源案例 × PNG/PDF/SVG：27 个文件通过，记录文件哈希与单页/有效 XML；包括三类坐标、三类标签、定量轨道及两个独立 gene/ribbon 层，变化 radius/curvature/gap/orientation。输入方向诊断为已有真实倒向 HSP，不是输出设备警告。
- Unicode 圆图 PDF 使用可工作的 Quartz，其他 PDF 使用标准设备；本机 Cairo 的 XQuartz 动态库不可用，未假称该设备通过。
- 文档：27 对入口、93 Rd、130 exports、865 段汉语参考；英语/汉语代码与图片共用；缺译或摘要过期阻止构建。21 张图片具有可复现源和哈希。
- 质粒图从旧会话记录对照专用主题/内部优先标签后，确认外侧混合布局仅在 >=30 酶切标签启用的缺口；现少量酶切和引物共存也按锚点共同排序。检查完整注释、引物与三个位点的导出图；非混合稀疏图路径保留。
- 按序列/链参数能力保留，独立双语教程和参数表；目录审计见 REPOSITORY-AUDIT-v0.16.5.md。主页 320 px，窄屏 220 px；README/长教程 200 px，无损原画不改。

本记录的本地通过不冒充远端部署。发布资产、正式状态和远端工作流最终以 GitHub Releases/Actions 为准；上线成功后再移除冻结 docs。源码包提供 SHA-256 校验和。v0.16.0–v0.16.4 是本次累积范围，未单独发布。


远端首轮 R-devel 的 162 测试通过，默认帮助 PDF 设备在 α 处失败；已为帮助
示例明确创建 ASCII 显示副本，并拆分长示例行。此修订不改运行时表达式或原始数据，
正式发布等待修订候选的远端检查结果。

源码包 SHA-256：`df099086bdde877a7e5c5ab87697bac909e8c0aa659fd70230a0716b372f9729`。
附件已下载回读并通过 SHA256SUMS.txt；发布时刻为 2026-11-02 00:27:49（Asia/Shanghai）。

## 远端与发布核对

- [四平台检查](https://github.com/DangJem/ggchord/actions/runs/41325331082)：Windows release、macOS release、Ubuntu release/devel 全部成功。
- [最终网站部署](https://github.com/DangJem/ggchord/actions/runs/41325563287)：源/产物检查、构建与部署成功。线上浏览器确认界面译文、语言切换、语言内帮助链接与放大 logo。
- R 包远端候选为 4d67400；最终站点提交 ba7f4c5 不改 R/man/DESCRIPTION/NAMESPACE/NEWS/vignette/data/inst。最终打包从 ba7f4c5 执行，本地标准 R CMD check 为 OK。
- [首次正式 Release](https://github.com/DangJem/ggchord/releases/tag/v0.16.5)：非草稿、非预发布，首次 tag 指向 ba7f4c5；附件服务端摘要与回读结果一致。

## README 修订与同版本重发（2026-11-02）

用户明确要求移除 README 中 reallinear.png 并重发 v0.16.5。两份 README 和生成源
同步删除指定图及专属图注，保留其他教程使用的图片资产；保留上一轮顶部 320 px
图标、逐图代码/结果布局及网站指南链接。标签按此次授权更新至 README 修订提交，
双语 Release 增补同等重发说明。首次发布事实不覆盖。

对比原标签，R、man、DESCRIPTION、NAMESPACE、NEWS、data、inst、vignettes 和
.Rbuildignore 均无差异。README、工具和版本记录不进入 R 源码包；因此继续使用
原来已验收且回读校验的源码包及 SHA256SUMS.txt，SHA-256 保持上述值。
本轮运行生成同步检查、双语文档源检查与 Git 差异检查，未重复运行完整发布矩阵。
