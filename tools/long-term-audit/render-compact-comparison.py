"""Read completed observations; no rendering or changes to source inputs."""
import argparse
import csv
import hashlib
import html
import json
from pathlib import Path

root = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--before', default='v0.27.0-cluster-final-results.json')
parser.add_argument('--after', default='v0.27.0-compact-final-results.json')
parser.add_argument('--output-prefix', default='v0.27.0-compact-comparison')
parser.add_argument('--image-retention')
args = parser.parse_args()
for name in (args.before, args.after, args.output_prefix):
    assert Path(name).name == name
before = json.loads((root / args.before).read_text())
after = json.loads((root / args.after).read_text())
retentions = [json.loads((root / 'v0.27.0-compact-image-retention.json').read_text())]
if args.image_retention:
    assert Path(args.image_retention).name == args.image_retention
    retentions.append(json.loads((root / args.image_retention).read_text()))
styles = args.before == 'v0.27.0-cluster-final-results.json' and args.after == 'v0.27.0-compact-final-results.json'
before_caption = Path(args.before).name.removesuffix('-results.json')
after_caption = Path(args.after).name.removesuffix('-results.json')


def image_file(file):
    name = Path(file).name
    reuse = next((r['reused_images'][name] for r in reversed(retentions) if name in r['reused_images']), None)
    if reuse:
        name = reuse['file']
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == reuse['sha256']
    assert (root / name).is_file(), name
    return html.escape(name, quote=True)


assert len(after['cases']) == 13, 'Do not publish a partial comparison as complete'
rows = []
sections = []
for name, case in after['cases'].items():
    old = before['cases'][name]
    for key in ('source_features', 'source_primers', 'selected_sites', 'selected_positions', 'site_policy'):
        assert case[key] == old[key], (name, key)
    assert [(r['source_row'], r['source_text'], r['component']) for r in case['visibility']] == [
        (r['source_row'], r['source_text'], r['component']) for r in old['visibility']]
    row = {'case': name}
    for key in ('n_visible_source_labels', 'n_hidden_source_labels', 'n_text_overlaps',
                'n_leader_crossings', 'n_leaders_through_text', 'n_label_outside_canvas',
                'mean_leader_length_mm', 'max_leader_length_mm'):
        row[key + '_before'] = old['quality'][key]
        row[key + '_after'] = case['quality'][key]
    row['converged_before'] = old['converged']
    row['converged_after'] = case['converged']
    row['anchors_off_own_feature_after'] = case['association']['anchors_off_own_feature']
    row['internal_glyphs_off_tangent_after'] = case['association']['internal_glyphs_off_tangent']
    row['inside_to_exterior_feature_labels'] = case['feature_exterior_from_inside']
    row['endpoint_strings'] = sum(a['endpoint'] for a in case['axis_labels'])
    row['hidden_endpoint_strings'] = sum(a['endpoint'] and not a['visible'] for a in case['axis_labels'])
    rows.append(row)
    desc = (f"可见来源 {case['quality']['n_visible_source_labels']}/{case['quality']['n_source_labels']}；"
            f"引线交叉 {old['quality']['n_leader_crossings']} → {case['quality']['n_leader_crossings']}；"
            f"穿字框 {old['quality']['n_leaders_through_text']} → {case['quality']['n_leaders_through_text']}；"
            f"收敛 {old['converged']} → {case['converged']}")
    sections.append(f'<section><h2>{html.escape(name)}</h2><p>{html.escape(desc)}</p><div class="pair">'
                    f'<figure><figcaption>之前：{html.escape(before_caption)}</figcaption><img src="{image_file(old["file"])}"></figure>'
                    f'<figure><figcaption>之后：{html.escape(after_caption)}，默认 individual</figcaption><img src="{image_file(case["file"])}"></figure>'
                    '</div></section>')
for name in (('pTRE_Tight_BI', 'pCAMBIA1300', 'pETDuet_1') if styles else ()):
    individual = image_file(f'v0.27.0-compact-{name}-individual.png')
    bracket = image_file(f'v0.27.0-compact-{name}-bracket.png')
    sections.append(f'<section><h2>{name}：显式样式对照</h2><p>本节不含限制位点；两侧保留相同完整特征与画布，只有 label_repeat 不同。</p><div class="pair">'
                    f'<figure><figcaption>individual</figcaption><img src="{individual}"></figure>'
                    f'<figure><figcaption>bracket</figcaption><img src="{bracket}"></figure></div></section>')
extra = root / 'v0.27.0-repeat-catalogue-control-results.json'
if styles and extra.exists():
    grouped = json.loads(extra.read_text())
    assert grouped['converged'] and all(row['visible'] for row in grouped['visibility'])
    sections.append('<section><h2>pTRIPZ：完整位点下的重复标注样式</h2>'
                    '<p>两侧相同原始特征、reference 位点及10×8英寸画布，6个连续 operator 共享一个计数标签，箭头保持独立。</p><div class="pair">'
                    f'<figure><figcaption>individual</figcaption><img src="{image_file("v0.27.0-compact-final-pTRIPZ.png")}"></figure>'
                    f'<figure><figcaption>bracket</figcaption><img src="{image_file("v0.27.0-repeat-pTRIPZ-bracket.png")}"></figure></div></section>')
with (root / (args.output_prefix + '.csv')).open('w', newline='') as stream:
    writer = csv.DictWriter(stream, rows[0].keys(), lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
(root / (args.output_prefix + '.html')).write_text('''<!doctype html><html lang="zh-Hans"><meta charset="utf-8">
<title>质粒紧凑布局：同输入、同尺寸对照</title><style>
body{font:16px system-ui,sans-serif;line-height:1.6;margin:2rem;background:#f5f5f5;color:#202020}
section{margin:2rem 0;padding:1rem;background:white;border-radius:8px}.pair{display:grid;grid-template-columns:1fr 1fr;gap:1rem}
figure{margin:0}img{width:100%;height:auto}figcaption{font-weight:bold}@media(max-width:800px){.pair{grid-template-columns:1fr}}
</style><h1>质粒紧凑布局：同输入、同尺寸对照</h1>
<p>前后均为10×8英寸、110 dpi，源DNA/特征/独立引物/位点profile相同。主对照均为individual默认样式。
新括注样式另列，不把共享来源计作独立字框。几何穿字框记录不等于最终可见墨迹穿越，也不等于已通过遮盖认证。
灰色虚线用于圆内特征→圈外标签；限制位点保留原样式。原点/终点文字优先保留；拥挤的普通中间刻度可依据重叠策略隐藏文字并保留刻度，诊断记录 hidden_near_backbone。源标签统计不包括坐标轴文字。</p>
''' + ''.join(sections) + '</html>')
print(f'Wrote {len(rows)} unchanged-input pairs and CSV; verified reused image hashes')
