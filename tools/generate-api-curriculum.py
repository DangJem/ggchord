#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Render the explicit API curriculum decisions; do not infer teaching depth."""
import argparse,html,json,re
from documentation_lib import ROOT,WEB

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');args=ap.parse_args()
    decisions=json.loads((WEB/'teaching-reviews.json').read_text())
    entries=json.loads((ROOT/'tools/teaching-audit/api-signatures.json').read_text())
    aliases={alias:p.stem for p in (ROOT/'man').glob('*.Rd') for alias in re.findall(r'\\alias\{([^}]+)\}',p.read_text())}
    assert set(decisions)=={e['name'] for e in entries},'Every actual export needs an explicit decision'
    for lang in ['en','zh']:
        zh=lang=='zh'
        lines=[('按 API 查找用途、主讲课程、模式、综合案例和完整参数。等级是教学范围的判断，不是所有参数组合的质量保证。L2 有可运行的最小使用与输入/输出解释；L3 还含关键参数与效果对照。上游重导出及纯拼写别名按家族说明，不重复制造课程。' if zh else 'Find purpose, primary lesson, modes, complete cases and full parameters by API. Levels assess teaching scope, not every parameter combination. L2 includes runnable minimal use and input/output explanation; L3 adds key-parameter contrasts. Upstream reexports and spelling aliases are explained by family rather than duplicate lessons.'),
          '\n## '+('API 学习路径' if zh else 'API learning paths')+' {#apis}\n',
          '| '+('API 与参考 | 用途 | 主讲与等级 | 模式/输入边界 | 参数图谱 | 综合案例/工作流' if zh else 'API/reference | Purpose | Lesson and level | Mode/input boundary | Parameter atlas | Case/workflow')+' |',
          '| --- | --- | --- | --- | --- | --- |']
        for e in entries:
            name=e['name'];r=decisions[name];path=WEB/lang/(r['page']+'.md')
            assert r['anchor'] in re.findall(r'\{#([^}]+)\}',path.read_text()),name+': missing anchor'
            purpose=html.escape(r['purpose'][lang]).replace('|','\\|')
            alias=('；同义拼写 '+r['alias_of'] if zh else '; spelling alias of '+r['alias_of']) if r.get('alias_of') else ''
            scope=('上游' if zh else 'upstream') if e['owner']!='ggchord' else r['level']
            lesson=f"{r['page']}.html#{r['anchor']}"
            lines.append(f"| [`{name}`](reference/{aliases[name]}.html) | {purpose}{alias} | [{scope}]({lesson}) | {r['modes_localized'][lang].replace('|','/')} | [↗](parameter-atlas.html#api-{name}) | [↗]({r['case']}.html) |")
        lines+=['\n## '+('怎样复现' if zh else 'How to reproduce')+' {#reproduce}\n',
         ('按主讲页从准备代码开始运行；字段含义在[基础表](chord-data.html)，实际格式转换在[导入](formats.html)。所有主讲页的共享片段在独立环境执行，重要对照图来自真实 ggchord 渲染。已标明的简化几何示意不作为生物学证据。完整版本、源码和验证记录随仓库保留；[安装页](installation.html)区分开发版与正式版。' if zh else 'Start at the preparation code in the primary lesson. Learn fields in [basic tables](chord-data.html) and conversions in [formats](formats.html). Shared lesson snippets execute in independent environments and important comparisons are rendered by ggchord. Explicit synthetic geometry demonstrations are not biological evidence. Keep version/source/verification records and use the [installation page](installation.html) to distinguish development and stable versions.')]
        result='\n'.join(lines)+'\n';path=WEB/lang/'api-curriculum.md'
        if args.check:assert path.read_text()==result,'Stale API curriculum: '+lang
        else:path.write_text(result)
if __name__=='__main__':main()
