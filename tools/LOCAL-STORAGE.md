# 本地生成物与空间维护

包运行使用源码与随包数据。旧安装副本、检查目录、网站预览和临时工具依赖不应
长期堆在 checkout；`.Rbuildignore` 已将这些开发目录排除出发布包。

```sh
# 只盘点，不删除
python3 tools/clean-local-artifacts.py

# 明确清理已核实的忽略生成物
python3 tools/clean-local-artifacts.py --clean
```

工具仅处理本地 `.site`、带检查日志的根级 `*.Rcheck`、旧 OCR 临时依赖，以及
`.codex-tmp/release-v*` / `cran-v*` 中已知的构建、安装、下载和设备输出子目录。
目标必须被 Git 忽略、没有跟踪文件且不是符号链接。清理前后核对所有原跟踪文件
内容与 Git refs；保留精简检查日志和网站来源快照到 `.codex-tmp/cleanup-records`。
恢复本地网站预览需显式重建 `.site`，清理不会部署或删除远端网站。

以下内容保留：`R/`、`data/`、原始 `examples/`、`man/` 及正式引用图像、冻结
`docs/`、设计文档、审计脚本与既有验证材料、个人备份和 Git 分支/标签/历史。
工具不清理整个 `.codex-tmp`，不执行 `git clean`、历史重写或自动垃圾回收。

日常诊断图、失败候选和大布局缓存应输出到系统临时目录或 `.codex-tmp`；需要
保留证据时，优先记录精简指标、来源和生成命令，只提交必要最终图。相同图像
按校验和复用已有文件。文档图仍在原位更新，避免为每次重跑另存整套图片。
旧发布日志、输入和来源记录不应为节省少量空间而一并清空。

本次实际清理与校验见
[v0.27.0-storage-cleanup-results.json](long-term-audit/v0.27.0-storage-cleanup-results.json)。
该记录是清理事实，不是包级发布认证。
