#!/usr/bin/env python3
"""Fresh-process lesson and version probes; no package build/check or testthat."""
import argparse,hashlib,json,re,subprocess
from documentation_lib import ROOT,WEB,snippets

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stable-library',required=True);ap.add_argument('--stable-archive',required=True);args=ap.parse_args()
    chunks=snippets();records=[]
    for label in ('stable','development'):
        setup=('library(ggchord, lib.loc='+json.dumps(args.stable_library)+')' if label=='stable' else 'pkgload::load_all(quiet=TRUE)')
        for slug in ('chord-start','collinear-start','plasmid-start'):
            text=(WEB/'en'/(slug+'.md')).read_text()
            keys=[key for group in re.findall(r'\{\{(?:code|prepare):([a-z,]+)\}\}',text) for key in group.split(',')]
            code=setup+'\n'+'\n'.join(chunks[k] for k in keys)
            objects=re.findall(r'^(p_\w+)\s*<-',code,re.M)
            code+='\ninvisible(ggplot2::ggplotGrob('+objects[-1]+'))\ncat(as.character(packageVersion("ggchord")))\n'
            version=subprocess.check_output(['Rscript','--vanilla','-e',code],cwd=ROOT,text=True).strip()
            records.append(dict(version=version,scope=label,page=slug,result='fresh R process and grob render succeeded'))
        probe=setup+'''\ncat(jsonlite::toJSON(list(
          version=as.character(packageVersion("ggchord")),
          nonadjacent_linewidth="nonadjacent_linewidth" %in% names(formals(geom_link_ribbon)),
          explicit_sequence_unit="sequence_unit" %in% names(formals(read_fasta_lengths)),
          protein_data=exists("protein_seq_data_example", inherits=TRUE)),auto_unbox=TRUE))
        p <- ggchord(plasmid_example_pUC19) + geom_seq() + coord_circular(rotation=0)
        file <- view_ggchord(p,width=6,height=6,units="in",viewer="none")
        file.copy(file,paste0("/tmp/ggchord-",'''+json.dumps(label)+''',"-circle.png"),overwrite=TRUE)
        '''
        output=subprocess.check_output(['Rscript','--vanilla','-e',probe],cwd=ROOT,text=True)
        first=output.index('{');data,end=json.JSONDecoder().raw_decode(output[first:]);records.append(dict(scope=label,probe=data,preview='/tmp/ggchord-'+label+'-circle.png'))
    result=dict(source_sha=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        stable_source='https://cran.r-project.org/src/contrib/ggchord_0.16.6.tar.gz',
        stable_archive_sha256=hashlib.sha256(__import__('pathlib').Path(args.stable_archive).read_bytes()).hexdigest(),records=records)
    (ROOT/'tools/teaching-audit/version-check.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print('PASS: three opening lessons rendered in separate processes on stable and development sources; version differences recorded')
if __name__=='__main__':main()
