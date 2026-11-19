#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the parameter inventory used by the handbook from current formals and English Rd."""
import argparse,importlib.util,json,re
from documentation_lib import ROOT,WEB
spec=importlib.util.spec_from_file_location('translations',ROOT/'tools/reference-translations.py')
ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
FAMILIES = ('coordinates', 'layers', 'positions', 'statistics', 'scales', 'guides', 'themes', 'workflows')


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

def kind(api, name, metadata):
    # Types are never inferred from a formal default or a shared spelling.
    entry = metadata.get(api, {}).get(name)
    if entry:
        return entry['type']['en']
    return 'Type, units and accepted forms: see this API argument contract at right'


def rd_arguments(aliases):
    contracts = {}
    for path in set(aliases.values()):
        text, fields = ref.inventory(path)
        values = {}
        for field in fields:
            if field['field'].startswith('argument:'):
                description = plain(text[field['start']:field['end']])
                for parameter in field['field'][9:].split(','):
                    values[parameter.strip()] = description
        contracts[path] = values
    return contracts


def argument_description(name, parameter, source):
    if name == 'geom_track_bar' and parameter == 'baseline':
        return ('Count bars require baseline 0; every nonzero value is rejected. '
                'Use geom_track_col() for supplied values and nonzero baselines.')
    return source.get(parameter, 'Defined by inherited arguments and upstream contract; see full reference.')


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
                    type=kind(name, parameter, overrides),
                    constraint=argument_description(name, parameter, contracts[path]),
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
                type=kind(name, parameter, overrides),
                constraint='Follow upstream help for types, units and accepted values.',
                source=str(path.relative_to(ROOT)) if path else entry['owner'],
                requiredness='upstream contract'
            )

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--check',action='store_true')
    ap.add_argument('--contracts-only',action='store_true',help='retained spelling for the parameter inventory (the default output)')
    args=ap.parse_args()
    entries=json.loads((ROOT/'tools/teaching-audit/api-signatures.json').read_text())
    aliases={a:p for p in (ROOT/'man').glob('*.Rd') for a in re.findall(r'\\alias\{([^}]+)\}',p.read_text())}
    overrides=json.loads((WEB/'content-en/parameter-overrides.json').read_text())
    contracts=rd_arguments(aliases)
    groups={fam:[entry for entry in entries if entry['owner']=='ggchord' and family(entry['name'])==fam] for fam in FAMILIES}
    upstream=[entry for entry in entries if entry['owner']!='ggchord']
    records=list(parameter_records(groups, upstream, aliases, contracts, overrides))
    target=ROOT/'tools/teaching-audit/parameter-contracts.json'
    result=json.dumps(records,ensure_ascii=False,indent=2)+'\n'
    if args.check:assert target.read_text()==result,'Stale parameter contracts'
    else:target.write_text(result)


if __name__ == '__main__':
    main()
