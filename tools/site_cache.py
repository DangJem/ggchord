#!/usr/bin/env python3
"""Content-addressed native intermediates; finalized pages are never cached."""
import argparse, hashlib, json, os, platform, re, shutil, subprocess, tempfile
from datetime import datetime, timezone
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def key(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def runtime():
    # Versions alone cannot detect replaced package data or native binaries.
    script = '''paths <- unique(vapply(.libPaths(),normalizePath,"")); files <- unlist(lapply(paths,function(p) list.files(p,recursive=TRUE,full.names=TRUE,all.files=TRUE))); files <- files[file.info(files)$isdir %in% FALSE]; hashes <- tools::md5sum(files); names(hashes) <- vapply(files,function(f) {i<-which(startsWith(f,paste0(paths,"/")))[1];paste0(i,":",substring(f,nchar(paths[i])+2))},""); cat(jsonlite::toJSON(list(R=unclass(R.version),home=R.home(),capabilities=capabilities(Xchk=FALSE),packages=as.list(hashes)),auto_unbox=TRUE))'''
    r = json.loads(subprocess.check_output(['Rscript','--vanilla','-e',script],cwd=ROOT,text=True))
    home=Path(r.pop('home'))
    engine={str(p.relative_to(home)):digest(p) for name in ['bin','lib','etc'] for p in (home/name).rglob('*') if p.is_file()}
    fonts = {}
    if shutil.which('fc-list') and platform.system() != 'Darwin':
        listing = subprocess.check_output(['fc-list','--format=%{file}\\n'],text=True)
        fonts = {name:digest(Path(name)) for name in sorted(set(listing.splitlines())) if Path(name).is_file()}
    elif platform.system() == 'Darwin':
        fonts = {str(p):digest(p) for base in ['/System/Library/Fonts','/Library/Fonts',str(Path.home()/'Library/Fonts')] for p in Path(base).rglob('*') if p.is_file()}
    else:
        raise RuntimeError('Font inventory unavailable: use --no-cache for a full build')
    native={}
    if shutil.which('ldconfig'):
        listing=subprocess.check_output(['ldconfig','-p'],text=True)
        native={name:digest(Path(name)) for name in sorted({line.split(' => ')[1] for line in listing.splitlines() if ' => ' in line}) if Path(name).is_file()}
    adapters={}
    for selector,default in [('ALIGNMENT_PYTHON','python3'),('DELTA_EXECUTABLE','show-coords')]:
        command=os.environ.get(selector,default); resolved=shutil.which(command)
        adapters[selector]=dict(command=command,file=digest(Path(resolved)) if resolved else None)
        if resolved and selector=='ALIGNMENT_PYTHON':
            probe="import sys,json,hashlib,importlib.metadata as m; print(json.dumps(dict(version=sys.version,files={str(d.metadata['Name'])+':'+str(f):hashlib.sha256(d.locate_file(f).read_bytes()).hexdigest() for d in m.distributions() for f in (d.files or []) if '..' not in f.parts and '__pycache__' not in f.parts and d.locate_file(f).is_file()}),sort_keys=True))"
            adapters[selector]['packages']=json.loads(subprocess.check_output([resolved,'-c',probe],text=True))
    return dict(R=r,engine=engine,fonts=fonts,native=native,adapters=adapters,os=platform.platform(),year=datetime.now(timezone.utc).year,
                pandoc=subprocess.check_output(['Rscript','--vanilla','-e','cat(as.character(rmarkdown::pandoc_version()))'],text=True),
                locale={'en':'en','zh':'zh_CN'},TZ='UTC',seed=1014,
                environment={k:v for k,v in os.environ.items() if k in {'R_LIBS_SITE','R_DEFAULT_PACKAGES','FONTCONFIG_FILE','FONTCONFIG_PATH','LC_ALL','LC_CTYPE','LANG'}})

def inputs(paths):
    return {str(p.relative_to(ROOT)):digest(p) for p in sorted(set(paths))
            if p.is_file() and '__pycache__' not in p.parts and
            not any(part.startswith('.') for part in p.relative_to(ROOT).parts)}

def shell_inputs():
    return inputs([ROOT/p for p in ['DESCRIPTION','NAMESPACE','LICENSE','LICENSE.md','_pkgdown.yml','tools/render-website.R','tools/build-website.py','tools/site_cache.py','website/pages.json','website/navigation.json']]+list((ROOT/'pkgdown').rglob('*')))

def reference_inputs():
    # Include new development inputs and installed teaching fixtures. Ignore
    # unrelated working copies and old images alongside repository examples.
    listed = subprocess.check_output(['git', 'ls-files', '-z', '--cached', '--others',
        '--exclude-standard', '--', 'examples', 'inst/extdata'], cwd=ROOT)
    fixtures = [ROOT/os.fsdecode(name) for name in listed.split(b'\0') if name]
    return inputs(fixtures + [p for name in ['R','man','data','website/reference-zh']
        for p in (ROOT/name).rglob('*')] + [ROOT/p for p in [
        'website/reference-examples.json','tools/reference-translations.py','tools/website-dependencies.R']])


def reference_cache_allowed():
    codes=list(json.loads((ROOT/'website/reference-examples.json').read_text()).values())
    for path in (ROOT/'man').glob('*.Rd'):
        text=path.read_text();at=text.find('\\examples{')
        if at>=0:
            start=at+len('\\examples{');i=start;depth=1
            while i<len(text) and depth:
                if text[i]=='\\':i+=2;continue
                depth+=(text[i]=='{')-(text[i]=='}');i+=1
            if depth:return False
            codes.append(text[start:i-1])
    for code in codes:
        if re.search(r'https?://|\b(?:Sys\.time|Sys\.Date|download\.file|socketConnection|url)\s*\(',code):return False
        selectors=re.findall(r'Sys\.getenv\(\s*"([^"\n]+)"',code)
        if len(selectors)!=len(re.findall(r'Sys\.getenv\s*\(',code)):return False
        if any(name not in {'ALIGNMENT_PYTHON','DELTA_EXECUTABLE'} for name in selectors):return False
    return True

class Cache:
    def __init__(self, out, enabled=True):
        self.out=out; self.enabled=enabled; self.root=ROOT/'.cache/website'; self.pending=ROOT/'.cache/pending-website'
        shutil.rmtree(self.pending,ignore_errors=True)
        self.report=[]; self.producer=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
        self.dirty=bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True))
    def restore(self,kind,identity,target):
        token=key(identity); entry=self.root/kind/token
        reason='disabled' if not self.enabled else 'not_found'
        if self.enabled and entry.is_dir() and not entry.is_symlink():
            try:
                receipt=json.loads((entry/'entry.json').read_text())
                assert receipt['schema']==1 and receipt['kind']==kind and receipt['identity']==identity and receipt['status']=='verified'
                assert not (entry/'files').is_symlink()
                assert receipt['files'] and all(not Path(n).is_absolute() and '..' not in Path(n).parts for n in receipt['files'])
                actual={str(p.relative_to(entry/'files')):digest(p) for p in (entry/'files').rglob('*') if p.is_file() and not p.is_symlink()}
                assert not any(p.is_symlink() for p in (entry/'files').rglob('*'))
                assert actual==receipt['files']
                target.mkdir(parents=True,exist_ok=True)
                shutil.copytree(entry/'files',target,dirs_exist_ok=True)
                self.report.append(dict(kind=kind,key=token,status='hit',producer=receipt['producer'],producer_dirty=receipt.get('dirty',False),language=identity.get('language'),page=identity.get('page',{}).get('slug',identity.get('slug'))))
                self.write(); return True
            except (AssertionError,KeyError,ValueError,OSError): reason='invalid_receipt_or_bytes'
        self.report.append(dict(kind=kind,key=token,status='miss',reason=reason,language=identity.get('language'),page=identity.get('page',{}).get('slug',identity.get('slug')))); self.write();return False
    def stage(self,kind,identity,source):
        if not self.enabled:return
        entry=self.pending/kind/key(identity);entry.mkdir(parents=True,exist_ok=True)
        assert not any(p.is_symlink() for p in source.rglob('*'))
        shutil.copytree(source,entry/'files',dirs_exist_ok=True)
        receipt=dict(schema=1,kind=kind,identity=identity,status='pending',producer=self.producer,dirty=self.dirty,
                     files={str(p.relative_to(source)):digest(p) for p in source.rglob('*') if p.is_file()})
        (entry/'entry.json').write_text(json.dumps(receipt,sort_keys=True)+'\n')
    def write(self):
        (self.out/'cache-report.json').write_text(json.dumps(self.report,indent=2)+'\n')

def publish():
    # Run the complete current site check before any entry becomes reusable.
    subprocess.run(['python3','tools/check-documentation.py'],cwd=ROOT,check=True)
    out=ROOT/'.site'; snap=json.loads((out/'source-snapshot.json').read_text())
    assert snap['status']=='verified unchanged across rendering and finalization'
    for group in ['render_inputs','archive_inputs']:
        assert all(digest(ROOT/n)==v for n,v in snap[group].items()),'Build inputs moved before cache promotion'
    manifest=json.loads((out/'build-manifest.json').read_text())
    for name, expected in manifest['source_sha256'].items():assert digest(ROOT/name)==expected
    for path in (ROOT/'.cache/pending-website').glob('*/*'):
        receipt=json.loads((path/'entry.json').read_text()); assert receipt['status']=='pending'
        assert {str(p.relative_to(path/'files')):digest(p) for p in (path/'files').rglob('*') if p.is_file()}==receipt['files']
        receipt['status']='verified';(path/'entry.json').write_text(json.dumps(receipt,sort_keys=True)+'\n')
        target=ROOT/'.cache/website'/receipt['kind']/key(receipt['identity']);target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists():shutil.rmtree(target)
        os.replace(path,target)
    shutil.rmtree(ROOT/'.cache/pending-website',ignore_errors=True)
    print('Verified current site and promoted native intermediate cache')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify-and-publish',action='store_true');a=ap.parse_args()
    if a.verify_and_publish:publish()
