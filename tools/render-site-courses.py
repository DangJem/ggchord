#!/usr/bin/env python3
"""Render only the explicitly registered new website PNGs, never a release matrix."""
import argparse,hashlib,json,subprocess,tempfile
from pathlib import Path
from documentation_lib import ROOT,WEB,snippets

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--only',nargs='+')
    args=parser.parse_args()
    rows=json.loads((WEB/'course-figures.json').read_text())
    if args.only:
        unknown=set(args.only)-{r['id'] for r in rows}
        if unknown:raise SystemExit('Unknown course figure IDs: '+str(unknown))
        rows=[r for r in rows if r['id'] in args.only]
    if not rows:raise SystemExit('No course figures selected')
    chunks=snippets()
    inputs=[*sorted((ROOT/'R').glob('*.R')),*sorted((ROOT/'data').glob('*.rda')),
            ROOT/'tools/documentation-examples.R',ROOT/'tools/site-course-examples.R',WEB/'course-figures.json']
    snapshot=lambda:{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    before=snapshot()
    jobs=[dict(row=r,code='\n\n'.join(chunks[k] for k in r['snippets'])) for r in rows]
    with tempfile.TemporaryDirectory(prefix='ggchord-courses-') as folder:
        folder=Path(folder);(folder/'jobs.json').write_text(json.dumps(jobs))
        script=folder/'render.R'
        script.write_text('''pkgload::load_all(quiet=TRUE)
args <- commandArgs(trailingOnly=TRUE)
jobs <- jsonlite::fromJSON(args[1], simplifyVector=FALSE)
for(job in jobs) {
  env <- new.env(parent=globalenv())
  eval(parse(text=job$code),env)
  plot <- get(job$row$object,env)
  ggplot2::ggsave(job$row$file, plot, width=job$row$width,
    height=job$row$height, units="in", dpi=job$row$dpi, bg="white")
  cat("RENDERED:",job$row$id,"\\n")
}
''')
        subprocess.run(['Rscript','--vanilla',str(script),str(folder/'jobs.json')],cwd=ROOT,check=True)
    if snapshot()!=before:raise SystemExit('Source changed while rendering; results not certified')
    path=ROOT/'tools/site-audit/course-rendering.json'
    previous=json.loads(path.read_text()) if path.exists() else {'figures':{}}
    for record in previous['figures'].values():
        record.setdefault('rendered_source_sha256',previous.get('source_sha256',{}))
    for job in jobs:
        row=job['row'];previous['figures'][row['id']]=dict(
            file=row['file'],width=row['width'],height=row['height'],dpi=row['dpi'],
            sha256=hashlib.sha256((ROOT/row['file']).read_bytes()).hexdigest(),
            code_sha256=hashlib.sha256(job['code'].encode()).hexdigest(),status='rendered PNG',
            rendered_source_sha256=before)
    previous.update(source_sha256=before,scope='Registered course PNGs at normal sizes; no package checks or release device matrix')
    path.write_text(json.dumps(previous,indent=2)+'\n')

if __name__=='__main__':main()
