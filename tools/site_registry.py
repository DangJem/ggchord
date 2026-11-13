"""Canonical page/section resolution shared by generation, search and migration checks."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
WEB=ROOT/'website'

def pages():return json.loads((WEB/'pages.json').read_text())

def canonical(target,context=None):
    if target.startswith(('../','http:','https:','mailto:')):return target
    if target.startswith('#'):
        if not context:return target
        target=context+'.html'+target
    if '.html' not in target:return target
    routes=json.loads((WEB/'page-migrations.json').read_text())
    visited=set()
    while target not in visited:
        visited.add(target)
        path,_,fragment=target.partition('#')
        mapping=routes.get(path,{})
        fresh=mapping.get(fragment)
        if fresh is None and not fragment:fresh=mapping.get('')
        if fresh is None or fresh==target:return target
        target=fresh
    raise ValueError('Migration cycle: '+str(visited))

def page_kind(path):
    match=next((p for p in pages() if p['canonical_path']==path and p['kind']!='alias'),None)
    return match['kind'] if match else 'api' if path.startswith('reference/') else 'utility'
