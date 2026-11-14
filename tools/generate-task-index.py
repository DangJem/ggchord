#!/usr/bin/env python3
"""Validate task-search destinations without recreating the removed overview page."""
import json,re
from documentation_lib import ROOT, WEB
from site_registry import canonical

def main():
    tasks=json.loads((WEB/'teaching-tasks.json').read_text())
    exports=set(re.findall(r'^export\(([^)]+)\)',(ROOT/'NAMESPACE').read_text(),re.M))
    for task in tasks:
        assert task['api'] in exports,task['api']
        old=task['page'];target=canonical(old.replace('#','.html#',1) if '#' in old else old+'.html')
        slug,_,fragment=target.partition('#');slug=slug.removesuffix('.html')
        for lang in ['en','zh']:
            source=WEB/lang/(slug+'.md');assert source.exists(),str(source)
            if fragment:
                text=source.read_text()
                assert '{#'+fragment+'}' in text or 'id="'+fragment+'"' in text,(old,target)
    assert not any((WEB/lang/'learn.md').exists() for lang in ['en','zh'])
    print('PASS: task-search records resolve to manual sections; removed overview stays absent')
if __name__=='__main__':main()
