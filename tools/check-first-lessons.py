#!/usr/bin/env python3
"""Run each first lesson in its own fresh R process, including a plot build."""
import re, subprocess
from documentation_lib import WEB, snippets
chunks=snippets()
for page in ('chord-start','collinear-start','plasmid-start'):
    text=(WEB/'en'/(page+'.md')).read_text()
    keys=[key for match in re.findall(r'\{\{(?:code|prepare):([\w,]+)\}\}',text) for key in match.split(',')]
    script='pkgload::load_all(quiet=TRUE)\n'+ '\n'.join(chunks[k] for k in keys)
    objects=re.findall(r'^(p_\w+)\s*<-',script,re.M)
    script+='\ninvisible(ggplot2::ggplot_build('+objects[-1]+'))\n'
    subprocess.run(['Rscript','--vanilla','-e',script],check=True,stdout=subprocess.DEVNULL)
    print('PASS fresh R process and plot build:',page)
