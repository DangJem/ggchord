#!/usr/bin/env python3
"""Check installed vignette URL closure in an explicitly built source tarball."""
import sys,tarfile,re,posixpath
from pathlib import Path
p=Path(sys.argv[1]);print('tarball bytes:',p.stat().st_size)
with tarfile.open(p) as t:
 names=set(t.getnames())
 for lang in ['ggchord_guide_hans','ggchord_vignette']:
  s=t.extractfile('ggchord/inst/doc/'+lang+'.html').read().decode();missing=[]
  for u in re.findall(r'(?:href|src)=["\x27]([^"\x27]+)',s):
   if u.startswith(('http:','https:','data:','#','mailto:','javascript:')):continue
   if posixpath.normpath(posixpath.join('ggchord/inst/doc',u.split('#')[0])) not in names:missing.append(u)
  print(lang,'missing:',missing);assert not missing
 print('PASS: installed offline resource closure')
