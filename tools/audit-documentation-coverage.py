#!/usr/bin/env python3
"""Inventory entry points and tutorial evidence; not semantic or runtime certification.

No package loading, example execution, vignette rendering or package checks.
"""
import csv
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

def main():
    entries = list(csv.DictReader((ROOT/'design-roadmap/audit-generated/candidate-exports-v0.16.4.csv').open()))
    # Parse actual calls, so prose mentions and commented calls are not examples.
    program = r'''
    expressions <- parse('tools/documentation-examples.R', keep.source=FALSE)
    calls <- character()
    walk <- function(x) {
      if (missing(x)) return(invisible(NULL))
      if (is.expression(x)) { for (e in x) walk(e); return(invisible(NULL)) }
      if (!is.call(x)) return(invisible(NULL))
      head <- x[[1L]]
      if (is.symbol(head)) calls <<- c(calls, as.character(head))
      if (is.call(head) && as.character(head[[1L]]) %in% c('::',':::'))
        calls <<- c(calls, as.character(head[[3L]]))
      for (e in as.list(x)[-1L]) walk(e)
      invisible(NULL)
    }
    walk(expressions)
    cat(jsonlite::toJSON(sort(unique(calls))))
    '''
    calls = set(json.loads(subprocess.check_output(['Rscript','-e',program],cwd=ROOT,text=True)))
    pages = {lang:{p.stem:p.read_text() for p in (ROOT/'website'/lang).glob('*.md')} for lang in ['en','zh']}
    output=[]
    for entry in entries:
        name = entry['name']
        pattern = r'(?<![\w.])'+re.escape(name)+r'(?![\w.])'
        mentions = {lang:[slug for slug,text in topics.items() if re.search(pattern,text)] for lang,topics in pages.items()}
        path = ROOT/entry['documentation']
        rd = path.read_text() if path.is_file() else ''
        output.append(dict(name=name,ownership=entry['ownership'],family=entry['family'],alias_of=entry['alias_of'],
            rd_source=entry['documentation'],rd_exists=path.is_file(),rd_has_examples='\\examples{' in rd,
            shared_code_has_call=name in calls,
            english_prose_pages=';'.join(sorted(mentions['en'])),han_prose_pages=';'.join(sorted(mentions['zh'])),
            interpretation='Reference presence, prose mentions and source calls are separate evidence; none certifies tutorial depth or execution'))
    path=ROOT/'design-roadmap/audit-generated/documentation-coverage-v0.16.5.csv'
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=output[0].keys(),lineterminator='\n');w.writeheader();w.writerows(output)
    owned=[row for row in output if row['ownership']=='ggchord']
    summary=dict(exports=len(output),owned=len(owned),upstream=len(output)-len(owned),
        reference_available=sum(row['rd_exists'] for row in output),
        owned_called_in_shared_source=sum(row['shared_code_has_call'] for row in owned),
        upstream_called_in_shared_source=sum(row['shared_code_has_call'] for row in output if row['ownership']!='ggchord'),
        scope='Source inventory only. Aliases and upstream objects do not each require an independent tutorial; task-level review remains necessary.')
    families=json.loads((ROOT/'website/coverage.json').read_text())
    summary['task_families']=len(families)
    summary['paired_entries']=len(pages['en'])
    summary['shared_figures']=len(json.loads((ROOT/'website/assets/captions.json').read_text()))
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
