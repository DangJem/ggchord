"""Render independently runnable excerpts from the shared R figure sources."""
from collections import defaultdict
from functools import lru_cache
from hashlib import sha256
import json
import re
import subprocess

from . import ROOT

RAW = re.compile(r'^```r\s*\n(.*?)^```[ \t]*$', re.M | re.S)
MACRO = re.compile(r'\{\{(code|prepare):([\w,]+)\}\}')
BLOCK = re.compile(RAW.pattern + '|' + MACRO.pattern, re.M | re.S)


def raw_key(code):
    return 'raw-' + sha256(code.strip().encode()).hexdigest()


@lru_cache(maxsize=1)
def example_index():
    from documentation_lib import snippets
    source = snippets()
    for page in sorted((ROOT / 'website/zh').glob('*.md')):
        for code in RAW.findall(page.read_text()):
            source[raw_key(code)] = code.strip()
    process = subprocess.run(
        ['Rscript', '--vanilla', str(ROOT / 'tools/website/example-index.R')],
        input=json.dumps(source), text=True, capture_output=True, check=True)
    parsed = json.loads(process.stdout)
    records = parsed['chunks']
    owners = defaultdict(list)
    for key, expressions in records.items():
        for expr in expressions:
            for symbol in expr['defines']:
                if key not in owners[symbol]:
                    owners[symbol].append(key)
    return source, records, owners, parsed['packages']


def independent_examples(text, localize, preview):
    source, records, owners, packages = example_index()
    page_keys = []

    def render(match):
        raw, macro, names = match.groups()
        keys = names.split(',') if names else [raw_key(raw)]
        if macro == 'prepare':
            page_keys.extend(keys)
            return ''
        # Sources inserted by a caller after reading the Han page (e.g. a
        # generated API reference) retain their explicitly authored R example.
        if any(k not in records for k in keys):
            return match[0]
        prior = list(page_keys)
        page_keys.extend(keys)
        active = set()
        emitted = set()
        pieces = []
        used_calls = set()
        plot_objects = set()
        last_is_plot = False

        def ensure(name, context, before):
            possible = owners.get(name, ())
            if not possible:
                return
            # Prefer the most recent definition in this example, then the page,
            # then a named shared source. Raw fragments never leak across pages.
            scopes = [context, *reversed(prior), *[k for k in possible if not k.startswith('raw-')]]
            for key in dict.fromkeys(scopes):
                if key not in possible:
                    continue
                cap = before if key == context else len(records[key])
                found = [i for i, item in enumerate(records[key][:cap]) if name in item['defines']]
                if not found:
                    continue
                for i in found:
                    emit(key, i)
                return

        def emit(key, i):
            nonlocal last_is_plot
            pair = (key, i)
            if pair in emitted or pair in active:
                return
            active.add(pair)
            record = records[key][i]
            for name in record['needs']:
                ensure(name, key, i)
            used_calls.update(record['calls'])
            lines = source[key].splitlines()
            part = '\n'.join(lines[record['start'] - 1:record['end']])
            last_is_plot = record['output'] in plot_objects or record['output'] in {
                'ggchord', 'ggplot', 'ggchord_fusion', 'fusion_magnify_auto', 'fusion_magnify_manual',
                'make_flagship_fusion', 'view_ggchord', 'ggsave'
            }
            if record['assigned']:
                for name in record['defines']:
                    if last_is_plot:
                        plot_objects.add(name)
                    else:
                        plot_objects.discard(name)
            if not part.startswith('library('):
                pieces.append(part)
            active.remove(pair)
            emitted.add(pair)

        for key in keys:
            for i in range(len(records[key])):
                emit(key, i)
            prior.append(key)
        if not pieces:
            return ''
        packages_needed = ['ggchord']
        # ggchord re-exports a subset of ggplot2. dplyr's data verbs (not
        # unrelated purrr/tidyr re-exports) take precedence for shared names.
        for pkg in ['ggplot2', 'dplyr', 'purrr', 'tidyr']:
            uses = used_calls.intersection(packages[pkg])
            if pkg == 'ggplot2':
                exported = re.findall(r'export\((\w+)\)', (ROOT / 'NAMESPACE').read_text())
                uses.difference_update(exported)
            if uses:
                packages_needed.append(pkg)
            used_calls.difference_update(packages[pkg])
        code = '\n'.join('library(' + p + ')' for p in packages_needed) + '\n\n'
        last = records[keys[-1]][-1]
        if last_is_plot and last['assigned'] and len(last['defines']) == 1:
            pieces.append(last['defines'][0])
        body = preview('\n'.join(pieces))
        code += localize(body)
        kind = 'plot' if last_is_plot or body.splitlines()[-1].startswith('view_ggchord(') else 'data'
        return '\n```{.r .run-example data-example-kind="' + kind + '"}\n' + code + '\n```\n'

    return BLOCK.sub(render, text)
