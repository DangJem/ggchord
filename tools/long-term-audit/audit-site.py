"""Static evidence only: page rows are review leads, not browser findings."""
import csv
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
PAGES = json.loads((ROOT / 'website/pages.json').read_text())


def write(name, rows, fields):
    with (OUT / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def prose_lines(text):
    fenced = False
    for number, line in enumerate(text.splitlines(), 1):
        if re.match(r'^\s*(```|~~~)', line):
            fenced = not fenced
            continue
        if not fenced:
            yield number, line


rows = []
issues = []
patterns = {
    'developer_status': r'本网站面向开发版|已实测的版本差异|尚未推送|本地.*(?:checkout|开发)|不要假定远端 main|源码 SHA|local (?:source )?checkout|not yet pushed|tested version differences|development version|local.*development|source SHA',
    'internal_sha': r'\b[0-9a-f]{40}\b',
    'hard_break_review': r'<br\s*/?>| {2}$',
}
for page in PAGES:
    for lang in ('zh', 'en'):
        path = ROOT / 'website' / lang / (page['slug'] + '.md')
        text = path.read_text() if path.exists() else ''
        found = []
        lines = list(prose_lines(text))
        for number, line in lines:
            for kind, pattern in patterns.items():
                if re.search(pattern, line, re.I):
                    found.append(kind)
                    issues.append(dict(path=str(path.relative_to(ROOT)), line=number,
                        kind=kind, evidence=line.strip(), interpretation='needs contextual review'))
        # 连续汉语源行仅提示检查，不能推断 HTML 强制换行。
        soft = []
        if lang == 'zh':
            for (n, a), (m, b) in zip(lines, lines[1:]):
                if m == n + 1 and re.search(r'[\u4e00-\u9fff]', a) and re.search(r'[\u4e00-\u9fff]', b) and not re.match(r'^\s*(?:[#|<>*-]|\d+[.])', a + b) and a.strip() and b.strip():
                    soft.append(n)
        rows.append(dict(slug=page['slug'], language=lang, title=page['title'][lang],
            group=page['group'], exists=path.exists(),
            source_hash=hashlib.sha256(text.encode()).hexdigest(),
            developer_status_hits=found.count('developer_status'),
            hard_break_review_hits=found.count('hard_break_review'),
            soft_line_candidates=';'.join(map(str, soft)),
            teaching_review='pending: purpose/prerequisites/input/reading/errors',
            rendering='not verified', browser='not verified: 320/375/650/768/1024/1280/1440',
            next_action='review paired content; render then inspect browser'))
write('website-pages.csv', rows, list(rows[0]))
write('website-findings.csv', issues, ['path','line','kind','evidence','interpretation'])

# 列出所有数值柱旧入口，不扫描冻结 docs 或历史设计。
calls = []
for folder in ('R', 'man', 'website', 'vignettes', 'tools', 'tests', 'examples'):
    for path in sorted((ROOT / folder).rglob('*')):
        if not path.is_file() or path.suffix not in ('.R', '.Rd', '.md', '.Rmd', '.py') or OUT in path.parents:
            continue
        for n, line in enumerate(path.read_text(errors='replace').splitlines(), 1):
            if 'geom_track_bar(' in line:
                calls.append(dict(path=str(path.relative_to(ROOT)), line=n, code=line.strip(),
                    migration='review numeric values -> col/identity; do not rename count examples blindly'))
write('track-bar-migration.csv', calls, ['path','line','code','migration'])

# 测试命中只是静态索引，不声称执行过。
exports = list(csv.DictReader((OUT / 'api.csv').open()))
references = []
for api in exports:
    for path in sorted((ROOT / 'tests/testthat').glob('*.R')):
        for n, line in enumerate(path.read_text().splitlines(), 1):
            if re.search(r'(?<![\w.])' + re.escape(api['name']) + r'\s*\(', line):
                references.append(dict(name=api['name'], path=str(path.relative_to(ROOT)), line=n,
                    evidence=line.strip(), execution='not run'))
write('test-references.csv', references, ['name','path','line','evidence','execution'])
print(f'{len(PAGES)} page pairs; {len(issues)} review leads; {len(calls)} bar references')

# 历史兼容取样只是证据索引，本轮未重新执行这些案例。
historical = list(csv.DictReader((ROOT / 'tools/teaching-audit/layer-compatibility.csv').open()))
capabilities = []
for api in exports:
    name = api['name']
    family = name.split('_')[0]
    if family not in ('geom', 'stat', 'position', 'coord'):
        continue
    responsibilities = dict(geom='几何与固定样式', stat='统计变换', position='局部放置', coord='全局坐标投影')
    role = ('track' if 'track' in name else 'annotation' if any(x in name for x in ('gene','feature','primer','restriction')) else
        'relationship' if any(x in name for x in ('link','ribbon')) else 'sequence')
    decision = '应扩展现有图层' if family == 'geom' else '可通过现有参数实现'
    if family == 'stat':
        decision = '需要 Stat'
    history = ';'.join(f"{r['coordinate']}:{r['status']}" for r in historical if r['api'] == name)
    capabilities.append(dict(capability=name, classification=decision, status='基线已有',
        responsibility=responsibilities[family], data_role=role, public_entry='保留现有入口',
        evidence=f"{api['path']}:{api['line']}", compatibility_evidence=history or '见函数契约；本轮未穷举',
        validation='静态源码；历史取样不等于本轮验收'))
candidates = [
 ('feature_data统一输入','只需新数据适配器','annotation','构造边界桥接；不新增生物学专用Geom'),
 ('SNP/单点/断点','可通过现有参数实现','annotation','marker/lollipop已有；VCF另作适配'),
 ('repeat/motif/domain/小InDel','应扩展现有图层','annotation','feature形状/角色/区间，不复制几何'),
 ('大型SV成对断点','应扩展现有图层','relationship','feature+pair/link；专用快捷入口待真实用例'),
 ('sequence logo','暂不开发','annotation','优先可选对接，不是简单区间形状'),
 ('统一geom_track','应新增通用层','track','源空间统计与Position后投影；快捷层薄包装'),
 ('track count/bin/summary/smooth/density/ecdf','需要 Stat','track','逐批证明口径，不能在投影后统计'),
 ('track point/path/step/area/ribbon/errorbar/text','应扩展现有图层','track','受控局部几何矩阵，不盲目透传'),
 ('纵向Collinear','可通过现有参数实现','sequence','扩展现有Coord轴向选项；该参数尚未实现'),
 ('镜像Sankey','需要新 Coord','relationship','左右实例注册；默认只过滤自连接，保留有向输入'),
 ('基础Dotplot','需要新 Coord','relationship','双轴投影，独立通用alignment mark薄层，无需Stat'),
 ('Radial spokes','需要新 Coord','relationship','线性辐条/循环邻接，只筛选输入HSP'),
 ('重复序列与between路由','应扩展现有图层','relationship','seq+link入口，共享实例桥接，不增Coord'),
 ('equal显示长度','可通过现有参数实现','sequence','一个opt-in策略，公开名字待审'),
 ('首末坐标轴','应扩展现有图层','sequence','显式scale设置优先，扩展axis测量/求解'),
 ('端点渐变/尖头','应扩展现有图层','relationship','共享ribbon几何，用户映射和标度优先'),
 ('pident分段和比例guide','可通过现有参数实现','relationship','stepsn与guide组合，需验证真实切割和非等宽图例'),
 ('质粒标签统一布局','应扩展现有图层','annotation','共享测量/solver，不能只调目录图参数'),
 ('adaptive布局/guide','可通过现有参数实现','sequence','渲染上下文与对象策略，公开名字待审'),
 ('Fusion/区域放大','应新增通用层','relationship','静态组合层，不是Coord'),
 ('VCF/BEDPE/MAF/chain/BigWig','只需新数据适配器','source','成熟生态可选依赖，不虚构缺失字段'),
 ('短序列字符/比对列','暂不开发','source','候选研究，真实字符和gap映射必需'),
 ('交互/在线绘图/Studio','暂不开发','none','明确非目标'),
]
for name, decision, role, reason in candidates:
    capabilities.append(dict(capability=name, classification=decision, status='目标或候选，未接入公开API',
        responsibility=reason, data_role=role, public_entry='候选；需里程碑验证', evidence='见DESIGN.md及prototype.R',
        compatibility_evidence='见设计规格', validation='方案；原型不等于公开实现'))
write('capabilities.csv', capabilities, list(capabilities[0]))
