#!/usr/bin/env python3
"""Conservative plot dependencies for reviewed, independent figure families.

Generation-file hashes remain historical provenance. Dependency hashes establish
whether the current engine and executable preparation are still equivalent.
Independent preparations and Fusion-only runtime files have explicit boundaries.
"""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def fingerprints():
    source = (ROOT / 'tools/documentation-examples.R').read_text()
    markers = list(re.finditer(r'^# ---- ([a-z]+) ----\n', source, re.M))
    chunks = {m[1]: source[m.end():markers[i + 1].start()].removesuffix('\n')
              for i, m in enumerate(markers[:-1])}
    # These chunks create only their own plots/tables. No other plot consumes
    # their outputs. Shared data, themes, imports and render helpers stay included.
    comparison = next(k for k, v in chunks.items() if 'label_plots <-' in v)
    runtime = sorted((ROOT / 'R').glob('*.R'))
    def runtime_digest(paths):
        return '\n'.join(str(p.relative_to(ROOT)) + ':' +
                         hashlib.md5(p.read_bytes()).hexdigest() for p in paths)
    data_digest = '\n' + '\n'.join(str(p.relative_to(ROOT)) + ':' +
        hashlib.md5(p.read_bytes()).hexdigest()
        for p in sorted((ROOT / 'data').glob('*.rda')) + [ROOT / 'R/sysdata.rda'])
    # Fusion implementation is entered only by Fusion objects. Shared render,
    # guide and label helpers still invalidate every figure as before.
    full_engine = runtime_digest(runtime) + data_digest
    engine = runtime_digest([p for p in runtime if p.name not in
                             {'fusion.R', 'fusion-render.R'}]) + data_digest
    captions = json.loads((ROOT / 'website/assets/captions.json').read_text())
    # Opening lessons and the chord showcase write separate objects. Existing
    # plots still depend conservatively on the shared source and runtime.
    starters = {'firstcomparison', 'firstchord', 'firstlinear', 'firstcircle', 'firstfeatures'}
    exercises = {'practicegenes', 'practicechord', 'practicelinear', 'practiceplasmid', 'firstsavebasic',
                 'inspectlinks', 'inspectgenes', 'inspectorigin'}
    practice_figures = {
        'practicechord': {'practicegenes', 'practicechord'},
        'practicelinear': {'practicegenes', 'practicelinear'},
        'practiceplasmid': {'practiceplasmid'},
    }
    front_chunks = {name for name in chunks if name.startswith('front')}
    front_dependencies = {
        'frontbackbone': {'frontbackbone'},
        'frontribbon': {'frontbackbone', 'frontlinks', 'frontribbon'},
        'frontfeature': {'frontbackbone', 'frontlinks', 'frontribbon', 'frontgenes', 'frontfeature'},
        'frontrotation': {'frontbackbone', 'frontrotation'},
        'frontgap': {'frontgap'},
    }
    for name in ('frontradius', 'frontcurvature', 'frontorientation', 'frontseqlabels', 'frontorder'):
        front_dependencies[name] = {name}
    for name in ('frontnoavoid', 'frontsmooth'):
        front_dependencies[name] = {'frontbackbone', 'frontlinks', 'frontobstacle', name}
    front_dependencies['frontroute'] = front_dependencies['frontribbon'] - {'frontribbon'} | {'frontroute'}
    front_dependencies['frontgenelabels'] = front_dependencies['frontfeature'] | {'frontgenelabels'}
    front_dependencies['frontmarkers'] = {'frontbackbone', 'frontlinks', 'frontmarkers'}
    linear_base = {'frontlinearbase'}
    linear_link = linear_base | {'frontlinks', 'frontlinearlink'}
    front_dependencies.update({
        'frontlinearbase': linear_base,
        'frontlinearlink': linear_link,
        'frontlinearcds': linear_link | {'frontgenes', 'frontlinearcds'},
        'frontlineargap': {'frontlinks', 'frontlineargap'},
        'frontlinearscale': linear_link | {'frontlinearscale'},
        'frontlinearcross': linear_link | {'frontlinearcross'},
        'frontcrosswidth': linear_link | {'frontlinearcross', 'frontcrosswidth'},
        'frontcrosssame': linear_link | {'frontlinearcross', 'frontcrosssame'},
        'frontcirclestyle': {'frontcirclestyle', 'firstcircle'},
        'frontcirclerotation': {'frontcirclerotation', 'firstcircle', 'firstfeatures'},
    })
    registry = json.loads((ROOT / 'website/parameter-comparisons.json').read_text())
    declared = {row['id']: set(row['snippets']) for row in registry}
    # Fusion preparations create only fusion_* / p_fusion* objects. Existing
    # figures do not consume them; record this boundary after a full export.
    declared.update({
        'fusionview': {'fusiondata', 'fusionpanels', 'fusionlinear', 'fusionview'},
        'fusionfour': {'fusiondata', 'fusionpanels', 'fusionlinear', 'fusionfour'},
    })
    independent = {name for name in chunks if name.startswith(('compare', 'fusion'))}
    result = {}
    for key in captions:
        dependency_engine = full_engine if key in {'fusionview', 'fusionfour'} else engine
        if key in declared:
            payload = dependency_engine + '\n' + '\n'.join(name + ':\n' + chunks[name]
                for name in chunks if name in declared[key])
            result[key] = hashlib.md5(payload.encode()).hexdigest()
            continue
        excluded = exercises - practice_figures.get(key, set())
        excluded.update(independent)
        # These new chunks write only independent comparison objects.
        if key not in {'alignleft', 'aligncenter', 'alignright', 'alignanchor', 'alignauto', 'tracktile'}:
            excluded.update({'alignbase', 'aligncenter', 'alignright', 'alignanchor', 'alignauto', 'tracktile'})
        excluded.update(front_chunks - front_dependencies.get(key, set()))
        if key not in starters and key not in practice_figures:
            excluded.update(starters - front_dependencies.get(key, set()))
        if key != 'showcasechord':
            excluded.update({'chordshowcasegenes', 'showcasechord'})
        if key != 'showcasecollinear':
            excluded.add('showcasecollinear')
        if key not in {'labelauto', 'labelradial', 'labelarc'}:
            excluded.add(comparison)
        payload = dependency_engine + '\n' + '\n'.join(
            name + ':\n' + value for name, value in chunks.items() if name not in excluded)
        result[key] = hashlib.md5(payload.encode()).hexdigest()
    return {'code': hashlib.md5(source.encode()).hexdigest(), 'assets': result}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stamp-current', action='store_true')
    args = parser.parse_args()
    result = fingerprints()
    if args.stamp_current:
        path = ROOT / 'website/assets/plots.csv'
        with path.open() as handle:
            reader = csv.DictReader(handle)
            fields, rows = reader.fieldnames, list(reader)
        # This migration is permitted only immediately after a complete current
        # export. Never infer equivalence from old PNGs or silently refresh code hashes.
        for row in rows:
            assert row['code_md5'] == result['code'], 'Full current export required'
            assert hashlib.md5((ROOT / row['file']).read_bytes()).hexdigest() == row['md5']
            row['dependency_md5'] = result['assets'][row['id']]
        if 'dependency_md5' not in fields:
            fields.append('dependency_md5')
        with path.open('w') as handle:
            writer = csv.DictWriter(handle, fields, lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
        print(f'Stamped reviewed dependencies for {len(rows)} current PNG exports.')
    else:
        print(json.dumps(result))


if __name__ == '__main__':
    main()
