# v0.16.6 CRAN incoming修正候选

2026-11-02；原GitHub发布标签不改，本轮修改保留v0.16.6本地版本分支。

## 原反馈与修正

- Windows/Debian原始日志均为4 NOTEs：incoming词语/包体积、Rd行宽、安装教程相对图片链接、帮助例子耗时；完整日志下载核对，与用户邮件一致。
- DESCRIPTION用Linear/reusable替代Collinear/composable，不通过WORDLIST规避incoming。原术语并非确证拼写错误。
- 两份教程共享48个离线图片（约1.88 MB）与SVG，脚本/样式内嵌，图片用data-external；.install_extras显式安装tutorial-assets和本地CSS。资源源于当前PNG，网站/README原图保留；man/figures不再进入源码包。
- 长帮助行拆分，3个慢帮助例子使用有来源的小输入，不用dontrun/donttest隐藏成功路径。四个帮助例子单独运行0.042/0.393/0.192/0.448秒，不将本机时长当平台门槛。
- roxygen/Rd同步；4项汉语译文的正文/参数完全不变后同步源摘要。47张PNG实际重建后全部与发布版逐字节相同，生成记录已更新。Bioconductor两类可选文档输入实测PASS。

## 验证与边界

- 两份安装HTML全部本地src/href在源码包内有目标；没有绝对机器路径或缺失的../man/figures链接。独立工具tools/check-vignette-resources.py保留此检查。
- 首轮联网as-cran的incoming feasibility OK，CRAN packages.rds发生一次60秒超时；其余原四类NOTE消失。发现的新figures目录NOTE已通过明确的tutorial-assets资源目录修正；本机旧Tidy的工具NOTE通过使用现有新版Tidy 5.8修正。
- 最终as-cran复查Status OK；完整Suggests、帮助例子、包内tests、双语vignette重建、PDF/HTML手册均通过。该轮仅关闭重复远端incoming查询（_R_CHECK_CRAN_INCOMING_REMOTE_=false），使用/opt/homebrew/bin新版Tidy，没有跳过测试或Suggests。
- 39对源文档、865译文、130导出/93 Rd检查通过。测试代码、运行时函数实现与公开签名未改。
- 全量check通过的候选为4,834,130字节；最后仅同步自动生成教程中可选依赖成功状态说明，再次R CMD build重建双语教程得到4,834,097字节，并重新验证安装资源闭包。未因这一状态文案重复整套check。
- 没有Windows/Debian新远端结果，不宣称CRAN已接受。不发送邮件、不代投webform、不重发GitHub、不移动tag、不合并main。

最终本地重投候选SHA-256：`0a5aaa55b589e86fa865694a8d065fd66a2a8e638b47e93e9e021f968e7195a7`。
原始包47,072,217字节；修正包4,834,097字节，约减少89.7%。
