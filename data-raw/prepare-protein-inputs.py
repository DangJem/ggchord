import tarfile,hashlib,json
from pathlib import Path
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("archive", type=Path)
args = parser.parse_args()
path = args.archive
root = Path(__file__).resolve().parent.parent
import os
os.chdir(root)
rows=[]
src=Path('examples/format-inputs/emales-mavirus-blastp.tsv')
for i,line in enumerate(src.read_text().splitlines(),1):
 if not line or line.startswith('#'): continue
 fields=line.split('\t')
 if fields[0]=='BVI_023A_0009':rows.append((float(fields[10]),-float(fields[11]),i,line,fields))
rows.sort()
_,_,line_num,line,fields=rows[0]
ids=fields[:2]

with tarfile.open(path) as tar:
 data={}
 for member in ['emales/emales.faa','emales/mavirus.faa']:
  original=tar.extractfile(member).read()
  records=original.decode().split('>')
  for record in records:
   if record and record.split()[0] in ids:
    data[record.split()[0]]=('>'+record,member,hashlib.sha256(original).hexdigest())
 assert len(data)==2
 out=Path('examples/format-inputs/protein-pair.faa')
 out.write_text(''.join(data[k][0] for k in ids))
 hits=Path('examples/format-inputs/protein-pair-blastp.tsv')
 hits.write_text('# BLASTP 2.9.0+\n'+line+'\n')
 provenance=json.loads(Path('examples/format-inputs/sources.json').read_text())
 provenance=[x for x in provenance if x.get('file') not in (str(out),str(hits))]
 provenance.append({'file':str(out),'url':'https://raw.githubusercontent.com/thackl/gggenomes/master/data-raw/emales.tgz',
 'archive_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
 'members':[{'member':data[k][1],'member_sha256':data[k][2],'id':k} for k in ids],
 'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),
 'transformation':'Two complete FASTA records selected by the query and subject IDs of the strongest BVI_023A_0009 BLASTp hit; headers and residues retained. Lengths are measured from complete proteins.'})
 provenance.append({'file':str(hits),'source_file':str(src),'source_line':line_num,
 'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'sha256':hashlib.sha256(hits.read_bytes()).hexdigest(),
 'transformation':'One observed HSP selected by ascending evalue then descending bitscore for BVI_023A_0009; original 12 fields retained; BLASTP program comment copied from the source.'})
 Path('examples/format-inputs/sources.json').write_text(json.dumps(provenance,indent=2)+'\n')
