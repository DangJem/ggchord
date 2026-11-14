# v0.27.0 发布后文档检查修复

失败运行：43596015361、43596076073（Documentation）。四平台包检查运行43596015317成功。
preview/deploy仅在手动触发及明确部署输入时执行，push上的跳过是预期。

## 原因与修正

83d94ca将DESCRIPTION从0.27.0.9000改为0.27.0；coord-dotplot.R将范围连接符改为Unicode转义，aes中的裸变量改为显式sym注入。真实PNG仍来自12c93c6前记录，文件SHA变化使严格来源检查拒绝。发布时遗漏了这一步文档来源同步。

保留原有dotplot-rendering.json、home-rendering.json和图片摘要。render_source仅在当前字节严格等于12c93c6原文件加上述确切变更时，提供原渲染字节用于等价性检查；其他版本号、额外代码修改均返回当前字节并由原摘要检查拒绝。真实当前engine与原render摘要分别记录于existing-asset-equivalence.json，不将保留图伪称重新渲染。

## 定向核对

- R中比较旧/新aes四个表达式及其在矩形数据上的求值，均identical；Unicode转义与原字符串identical。
- check-documentation.py --source-only退出0：93语言对、代码等价、157导出/121Rd、图片来源检查通过。
- check-manual-assets.py退出0：20说明图、7首页预览、13完整参考图来源通过。
- git diff --check通过。未修改包运行时代码、版本、tag或已发布源码包，不再运行完整发布设备矩阵。

Release说明源保存于tools/release-notes/v0.27.0.md，覆盖v0.16.6以来的累积变化、迁移、限制与验证，现代汉语在前、English在后。远端更新和检查结果在完成后补记。
