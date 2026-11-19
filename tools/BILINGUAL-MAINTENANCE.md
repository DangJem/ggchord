# 双语内容维护约定

今后维护双语网站、README、指南、函数帮助和其他成对说明材料时，先写英语原文，
再翻译并校对现代汉语。修改英语后同步译文。共享代码、数据、图片、章节结构、
锚点和导航标识。双语 Release 仍以现代汉语在前、English 在后的顺序展示。
既有说明在下次维护时按此规则更新，历史发布记录保留当时的内容。
这份 Git 跟踪的约定适用于后续版本及新的项目 checkout。

以 `website/en/` 的英语 Markdown 和 `website/content-en/` 的结构化英语内容为源。
逐段校对汉语译文后，保存在 `website/translations/zh/`。页面标识、代码、图像、
链接目标和章节锚点共用。分别由两种语言的 Markdown 生成 HTML。

`website/content-policy.json` 按条目登记源、译文与输出。英语改动后执行
`python3 tools/sync-website-languages.py --prepare`，翻译并核对标为
`needs_review` 的字段，更新状态后运行同步器和 `--check`。整页及结构化对象
可分别使用 `--reviewed-page`、`--reviewed-content` 提交已校对稿。

README 的版式、链接和 R 代码维护在 `website/readme-template.md`，英语文案维护在
`website/content-en/readme-text.json`。生成配对 JSON 后，用
`python3 tools/sync-tutorials.py` 更新两份 README。独立离线教程的英语源和已校对
汉语稿在 `vignettes/` 下。完整操作示例和构建说明见[网站维护指南](../website/README.md)。
