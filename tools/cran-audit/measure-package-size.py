# Read-only source/data size comparison; never creates a package archive.
# Run from the source root; optional argument selects the prior source commit.
import json,subprocess,tarfile,hashlib,sys
from pathlib import Path
root_files={'DESCRIPTION','LICENSE','NAMESPACE','NEWS.md'}
def eligible(n):
 return n in root_files or (n.split('/')[0] in {'R','data','man','inst','tests','vignettes'} and not n.startswith(('man/figures/','inst/doc/')))
def manifest(ref):
 out={}
 for line in subprocess.check_output(['git','ls-tree','-rl',ref],text=True).splitlines():
  meta,n=line.split('\t',1)
  if eligible(n):out[n]=int(meta.split()[3])
 return out
baseline=manifest('refs/tags/v0.28.1')
prior_ref=sys.argv[1] if len(sys.argv)>1 else 'HEAD'
prior=manifest(prior_ref)
current={n:Path(n).stat().st_size for n in subprocess.check_output(['git','ls-files'],text=True).splitlines() if eligible(n) and Path(n).is_file()}
p=Path('ggchord_0.28.1.tar.gz')
with tarfile.open(p) as t:
 members={m.name.split('/',1)[1]:m.size for m in t if m.isfile() and '/' in m.name}
size=dict(baseline_tag='refs/tags/v0.28.1',baseline_commit=subprocess.check_output(['git','rev-parse','refs/tags/v0.28.1^{commit}'],text=True).strip(),prior_turn_commit=subprocess.check_output(['git','rev-parse',prior_ref],text=True).strip(),existing_archive_bytes=p.stat().st_size,existing_archive_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),existing_archive_payload_bytes=sum(members.values()),data=dict(baseline_bytes=sum(v for n,v in baseline.items() if n.startswith('data/')),current_bytes=sum(v for n,v in current.items() if n.startswith('data/')),added_files={n:v for n,v in current.items() if n.startswith('data/') and n not in baseline}),tracked_package_source_bytes=dict(baseline=sum(baseline.values()),prior_turn=sum(prior.values()),current=sum(current.values())),notes='No new archive built. Tracked package source scope: R,data,man (excluding figures),inst (excluding generated doc),tests,vignettes and DESCRIPTION/LICENSE/NAMESPACE/NEWS. Excludes dev-only website,tools,data-raw and frozen docs per .Rbuildignore. Generated inst/doc and build metadata require an actual authorized package build to determine compressed size.')
size['data']['delta_bytes']=size['data']['current_bytes']-size['data']['baseline_bytes']
size['tracked_package_source_bytes']['cumulative_delta']=sum(current.values())-sum(baseline.values())
size['tracked_package_source_bytes']['this_turn_delta']=sum(current.values())-sum(prior.values())
Path('tools/cran-audit/package-size-v0.28.1.json').write_text(json.dumps(size,indent=2)+'\n')
print(json.dumps(size,indent=2))
