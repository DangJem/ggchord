#!/usr/bin/env python3
"""Select unchanged public records; never execute or broadly extract an archive."""
import hashlib,json,tarfile,argparse
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('archive',type=Path);a=p.parse_args()
root=Path(__file__).resolve().parent.parent;out=root/'examples/format-inputs';out.mkdir(exist_ok=True)
archive_hash=hashlib.sha256(a.archive.read_bytes()).hexdigest()
metadata=json.loads((out/'sources.json').read_text())
metadata=[x for x in metadata if x.get('accession')=='L09137.2' or Path(x.get('file','')).name.startswith('protein-pair')]
with tarfile.open(a.archive) as tar:
 def select(member,name,rule):
  original=tar.extractfile(member).read();lines=original.decode().splitlines(keepends=True)
  selected=[(i,l) for i,l in enumerate(lines,1) if rule(l)]
  data=''.join(l for _,l in selected).encode();path=out/name;path.write_bytes(data)
  ranges=[]
  for i,_ in selected:
   if ranges and i==ranges[-1][1]+1:ranges[-1][1]=i
   else:ranges.append([i,i])
  metadata.append(dict(file=str(path.relative_to(root)),url='https://github.com/thackl/gggenomes/raw/master/data-raw/emales.tgz',archive_sha256=archive_hash,member=member,member_sha256=hashlib.sha256(original).hexdigest(),sha256=hashlib.sha256(data).hexdigest(),source_line_ranges=ranges,retrieved='2026-11-01',rule=rule.__doc__ or 'Unchanged complete member'))
 def paf(l):
  """Keep BVI_023A -> Cflag_017A records, including original optional tags."""
  f=l.split('\t');return len(f)>=12 and f[0]=='BVI_023A' and f[5]=='Cflag_017A'
 def gff(l):
  """Keep original GFF headers and BVI_023A annotations."""
  return l.startswith('#') or l.startswith('BVI_023A\t')
 select('emales/emales.paf','emales-pair.paf',paf)
 select('emales/emales_mavirus-blastp.tsv','emales-mavirus-blastp.tsv',lambda l:True)
 select('emales/emales.gff','emales-BVI_023A.gff3',gff)
 select('emales/emales.fna.split/emales.id_BVI_023A.fna','emales-BVI_023A.fna',lambda l:True)
(out/'sources.json').write_text(json.dumps(metadata,indent=2)+'\n')
print('Selected real PAF, BLASTp, GFF3 and FASTA records with hashes and physical source lines.')
