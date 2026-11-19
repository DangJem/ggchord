#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Current signatures plus canonical, already-translated Rd argument contracts."""
import argparse,html,importlib.util,json,re
from documentation_lib import ROOT,WEB
spec=importlib.util.spec_from_file_location('translations',ROOT/'tools/reference-translations.py')
ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
TITLES={
    'coordinates':('坐标与全局视图','Coordinates and global views'),
    'layers':('几何与标签','Geometry and labels'),
    'positions':('局部放置','Local placement'),
    'statistics':('区间和连接统计','Interval and link statistics'),
    'scales':('角色标度','Role scales'),
    'guides':('图例与 guide','Legends and guides'),
    'themes':('主题与物理外观','Themes and physical appearance'),
    'workflows':('输入、转换、诊断与输出','Input, conversion, diagnostics and output')
}

def plain(text):
    text=re.sub(r'\\link(?:\[[^]]*\])?\{([^{}]*)\}',r'\1',text)
    for _ in range(8):
        text=re.sub(r'\\[A-Za-z]+(?:\[[^]]*\])?\{([^{}]*)\}',r'\1',text)
    return re.sub(r'\s+',' ',text.replace('\\%','%').replace('\\_','_')).strip()

def family(name):
    if name == "layout_sequences":
        return "coordinates"
    if name == "track_band":
        return "positions"
    for prefix,key in [('scale_','scales'),('theme_','themes'),('guide_','guides'),('position_','positions'),('geom_','layers'),('coord_','coordinates'),('stat_','statistics')]:
        if name.startswith(prefix):return key
    return 'workflows'

def kind(api,name,lang,metadata):
    # Types are never inferred from a formal default or a shared spelling.
    # Uncurated entries defer to their exact canonical Rd argument contract.
    entry=metadata.get(api,{}).get(name)
    if entry:return entry['type'][lang]
    return ('类型、单位及允许格式见右侧本函数参数契约' if lang=='zh' else
            'Type, units and accepted forms: see this API argument contract at right')

def rd_arguments(aliases, translated):
    contracts={}
    for path in set(aliases.values()):
        text,fields=ref.inventory(path)
        saved=json.loads((WEB/'reference-zh'/(path.stem+'.json')).read_text())['translations'] if translated else None
        values={}
        for field in fields:
            if not field['field'].startswith('argument:'):continue
            descriptions={'en':plain(text[field['start']:field['end']])}
            if saved is not None:
                zh=saved[field['key']]
                for i,token in enumerate(field['tokens']):zh=zh.replace('⟦'+str(i)+'⟧',token)
                descriptions['zh']=plain(zh)
            for parameter in field['field'][9:].split(','):values[parameter.strip()]=descriptions
        contracts[path]=values
    return contracts

def argument_description(name, parameter, source, lang):
    if name == 'geom_track_bar' and parameter == 'baseline':
        if lang == 'zh':
            return '计数柱仅允许 0 基线；任何非零值均报错。给定数值和非零基线使用 geom_track_col()。'
        return ('Count bars require baseline 0; every nonzero value is rejected. '
                'Use geom_track_col() for supplied values and nonzero baselines.')
    fallback = ('由继承参数及上游接口定义；见完整参考。' if lang == 'zh' else
                'Defined by inherited arguments and upstream contract; see full reference.')
    return source.get(parameter, {}).get(lang, fallback)

def parameter_records(groups, upstream, aliases, contracts, overrides):
    for group in groups.values():
        for entry in group:
            name = entry['name']
            path = aliases[name]
            for parameter, default in (entry['parameters'] or {}).items():
                type_source = ('explicit API override' if parameter in overrides.get(name, {}) else
                               'canonical Rd argument contract')
                yield dict(
                    api=name, parameter=parameter, default=default,
                    type=kind(name, parameter, 'en', overrides),
                    constraint=argument_description(name, parameter, contracts[path], 'en'),
                    source=str(path.relative_to(ROOT)),
                    requiredness='see API contract; not inferred from formals',
                    type_source=type_source, units='see canonical argument contract'
                )
    for entry in upstream:
        name = entry['name']
        path = aliases.get(name)
        for parameter, default in (entry['parameters'] or {}).items():
            yield dict(
                api=name, parameter=parameter, default=default,
                type=kind(name, parameter, 'en', overrides),
                constraint='Follow upstream help for types, units and accepted values.',
                source=str(path.relative_to(ROOT)) if path else entry['owner'],
                requiredness='upstream contract'
            )

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--check',action='store_true')
    ap.add_argument('--contracts-only',action='store_true',help='refresh only the developer parameter contract JSON')
    args=ap.parse_args()
    entries=json.loads((ROOT/'tools/teaching-audit/api-signatures.json').read_text())
    aliases={a:p for p in (ROOT/'man').glob('*.Rd') for a in re.findall(r'\\alias\{([^}]+)\}',p.read_text())}
    overrides=json.loads((WEB/'parameter-overrides.json').read_text())
    contracts=rd_arguments(aliases, translated=not args.contracts_only)
    groups={fam:[entry for entry in entries if entry['owner']=='ggchord' and family(entry['name'])==fam] for fam in TITLES}
    upstream=[entry for entry in entries if entry['owner']!='ggchord']
    records=list(parameter_records(groups, upstream, aliases, contracts, overrides))
    if not args.contracts_only: write_atlas(groups, upstream, aliases, contracts, overrides, args.check)
    target=ROOT/'tools/teaching-audit/parameter-contracts.json'
    result=json.dumps(records,ensure_ascii=False,indent=2)+'\n'
    if args.check:assert target.read_text()==result,'Stale parameter contracts'
    else:target.write_text(result)
    if not args.contracts_only:
        import subprocess
        subprocess.run(['python3','tools/sync-website-languages.py',*(['--check'] if args.check else [])],cwd=ROOT,check=True)

def write_atlas(groups, upstream, aliases, contracts, overrides, check):
    from website.source_links import source_document

    migrations = json.loads((WEB/'page-migrations.json').read_text()).get('parameter-atlas.html', {})
    active_ids = {'api-'+entry['name'] for group in groups.values() for entry in group}
    active_ids.update('api-'+entry['name'] for entry in upstream)
    retired_anchors = {anchor: url for anchor, url in migrations.items()
                       if anchor.startswith('api-') and anchor not in active_ids}

    bodies = {}
    for lang in ('zh', 'en'):
        zh=lang=='zh'
        content=[('本页按用途汇总当前完整函数参数，不替代用法说明。默认表达式来自实际 formals，说明来自同版本 Rd；NULL 可能在布局阶段解析为默认值。候选向量既可能供 match.arg 选择，也可能表示同时启用多个选项，须阅读说明。上游重导出遵循其自身帮助。' if zh else 'Browse complete current parameters by purpose. Defaults come from actual formals and descriptions from matching Rd. NULL can resolve during layout. A choice vector may mean match.arg alternatives or multiple enabled options; read its description. Upstream reexports follow their own help.'),
         ('\n优先从[功能索引](reference/index.html)进入用法，重要参数对照见对应说明，完整签名和额外参数转发见函数参考。以下按需展开；不会把上游 ... 的全部可能参数伪装成 ggchord 的独立 API。\n' if zh else '\nStart with the [function index](reference/index.html). Examples compare important parameters; references define full signatures and forwarding. Expand functions as needed; upstream dots are not presented as separate ggchord APIs.\n')]
        content += ['<label class="atlas-filter">'+('查找函数或参数' if zh else 'Find a function or parameter')+'<input type="search" id="parameter-filter" autocomplete="off"></label><p id="parameter-filter-status" class="atlas-status" role="status"></p>\n']
        for fam,group in groups.items():
            content+=['\n## '+TITLES[fam][0 if zh else 1]+' {#'+fam+'}\n']
            for entry in group:
                name=entry['name']
                path=aliases[name];source=contracts[path]
                content+=['<span id="api-'+name+'"></span>\n<details><summary><code>'+name+'()</code></summary>\n',
                    ('\n[完整参考]' if zh else '\n[Full reference]')+'(reference/'+path.stem+'.html)\n',
                    ('\n| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |\n| --- | --- | --- | --- |' if zh else '\n| Parameter | Signature default | Type/value rule | Meaning, range and constraints |\n| --- | --- | --- | --- |')]
                for parameter,default in (entry['parameters'] or {}).items():
                    description=argument_description(name,parameter,source,lang)
                    default_text=default or ('无形式默认；是否必需见契约' if zh else 'No formal default; requiredness follows the contract')
                    cells=[f'`{parameter}`',f'`{default_text}`',kind(name,parameter,lang,overrides),description]
                    content+=['| '+' | '.join(c.replace('|','\\|').replace('\n',' ') for c in cells)+' |']
                content+=['\n</details>\n']
        content += ['\n## '+('上游重导出与对象' if zh else 'Upstream reexports and objects')+' {#upstream}\n',
                    ('这些入口保留上游 ggplot2/grid 语义；它们不是新增的 ggchord 图层。after_stat/after_scale 的阶段不同，不能由 identical 结果推断可互换。\n' if zh else 'These retain upstream ggplot2/grid semantics rather than defining new ggchord layers. After_stat and after_scale differ by stage even if identical compares their R objects equal.\n')]
        for entry in upstream:
            name=entry['name'];path=aliases.get(name)
            content += ['<span id="api-'+name+'"></span>\n<details><summary><code>'+name+'</code> · '+html.escape(entry['owner'])+'</summary>\n']
            if path:content += [('上游入口参考' if zh else 'Upstream entry reference')+': [Rd](reference/'+path.stem+'.html)\n']
            content += ['| '+('参数 | 签名默认 | 类型/约束' if zh else 'Parameter | Signature default | Type/contract')+' |','| --- | --- | --- |']
            for parameter,default in (entry['parameters'] or {}).items():
                note=('遵循上游接口；具体类型、单位和允许值见其帮助。' if zh else 'Follow upstream help for types, units and accepted values.')
                content += ['| `'+parameter+'` | `'+(default or ('无形式默认；见上游契约' if zh else 'No formal default; see upstream contract'))+'` | '+kind(name,parameter,lang,overrides)+'; '+note+' |']
            content += ['\n</details>\n']
        if retired_anchors:
            title = '旧参数地址的替代条目' if zh else 'Current entries for old parameter links'
            label = '查看当前用法' if zh else 'View current usage'
            content.extend(['<noscript>', '<nav aria-labelledby="retired-parameter-links">',
                            '<p id="retired-parameter-links">'+title+'</p>', '<ul>'])
            for anchor, url in retired_anchors.items():
                identity = html.escape(anchor, quote=True)
                name = html.escape(anchor.removeprefix('api-'))
                href = html.escape(url, quote=True)
                content.append(f'<li id="{identity}"><code>{name}()</code> · <a href="{href}">{label}</a></li>')
            content.extend(['</ul>', '</nav>', '</noscript>'])
        bodies[lang] = source_document('\n'.join(content).rstrip()+'\n')
    from website.generated_manual import write_page
    write_page('parameter-atlas', bodies, check)
if __name__=='__main__':main()
