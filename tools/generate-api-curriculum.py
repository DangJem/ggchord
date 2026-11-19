#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Render the explicit API curriculum decisions; do not infer teaching depth."""
import argparse,html,json,re
from documentation_lib import ROOT,WEB
from website.source_links import source_document

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');args=ap.parse_args()
    decisions=json.loads((WEB/'teaching-reviews.json').read_text())
    entries=json.loads((ROOT/'tools/teaching-audit/api-signatures.json').read_text())
    aliases={alias:p.stem for p in (ROOT/'man').glob('*.Rd') for alias in re.findall(r'\\alias\{([^}]+)\}',p.read_text())}
    assert set(decisions)=={e['name'] for e in entries},'Every actual export needs an explicit decision'
    bodies = {}
    for lang in ['zh', 'en']:
        zh=lang=='zh'
        lines=[('按函数名查找用途、适用输入和示例。函数名链接到签名和参数说明。' if zh else 'Look up functions by name, purpose and accepted input. Each function name links to its signature and argument descriptions.'),
          '\n## '+('按函数名查找' if zh else 'Find a function')+' {#apis}\n',
          '| '+('API 与参考 | 用途 | 用法说明 | 模式/输入边界 | 参数图谱 | 综合案例/工作流' if zh else 'API/reference | Purpose | Usage section | Mode/input boundary | Parameter atlas | Case/workflow')+' |',
          '| --- | --- | --- | --- | --- | --- |']
        for e in entries:
            name=e['name'];r=decisions[name];path=WEB/lang/(r['page']+'.md')
            assert r['anchor'] in re.findall(r'\{#([^}]+)\}',(WEB/'zh'/(r['page']+'.md')).read_text()),name+': missing anchor'
            purpose=html.escape(r['purpose'][lang]).replace('|','\\|')
            alias=('；同义拼写 '+r['alias_of'] if zh else '; spelling alias of '+r['alias_of']) if r.get('alias_of') else ''
            scope=('上游' if zh else 'upstream') if e['owner']!='ggchord' else ('用法' if zh else 'Usage')
            lesson=f"{r['page']}.html#{r['anchor']}"
            lines.append(f"| [`{name}`](reference/{aliases[name]}.html) | {purpose}{alias} | [{scope}]({lesson}) | {r['modes_localized'][lang].replace('|','/')} | [↗](parameter-atlas.html#api-{name}) | [↗]({r['case']}.html) |")
        lines += ['\n## '+('运行示例' if zh else 'Run an example')+' {#reproduce}\n',
                  ('打开对应条目，复制要绘制的图或要执行的操作。[查看数据格式](formats.html)。' if zh else
                   'Open the entry and copy the example for the figure or operation you need. [Browse data formats](formats.html).')]
        bodies[lang] = source_document('\n'.join(lines).rstrip()+'\n')
    from website.generated_manual import write_page
    write_page('api-curriculum', bodies, args.check)
    import subprocess
    subprocess.run(['python3','tools/sync-website-languages.py',*(['--check'] if args.check else [])],cwd=ROOT,check=True)
if __name__=='__main__':main()
