#!/usr/bin/env python3
"""Run each first lesson in its own fresh R process, including a plot build."""
import re, subprocess, json, hashlib
from pathlib import Path
from documentation_lib import WEB, snippets
chunks=snippets()
records=[]
for page in ('multi-sequence-start','circular-start','dotplot-start','fusion-start'):
    text=(WEB/'en'/(page+'.md')).read_text()
    keys=[key for match in re.findall(r'\{\{(?:code|prepare):([\w,]+)\}\}',text) for key in match.split(',')]
    script='pkgload::load_all(quiet=TRUE)\n'+ '\n'.join(chunks[k] for k in keys)
    objects=re.findall(r'^(p_\w+)\s*<-',script,re.M)
    script+='\ninvisible('+('export_ggchord_fusion('+objects[-1]+', width=12, height=8)' if page=='fusion-start' else 'ggplot2::ggplot_build('+objects[-1]+')')+')\n'
    wrapped='invisible(eval(parse(text='+json.dumps(script)+')))'
    subprocess.run(['Rscript','--vanilla','-e',wrapped],check=True)
    print('PASS fresh R process and plot build:',page)
    records.append(dict(page=page,result='fresh R process and public build/export succeeded',
        script_sha256=hashlib.sha256(script.encode()).hexdigest()))
Path('tools/site-audit/first-lessons.json').write_text(json.dumps(dict(
    records=records,scope='Four independent opening lessons; no install/network, package tests or device matrix'),indent=2)+'\n')
