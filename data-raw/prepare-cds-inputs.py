#!/usr/bin/env python3
"""Extract real CDS evidence without changing coordinates or source characters."""
from pathlib import Path
import hashlib
import json
from Bio import SeqIO
import Bio

root = Path(__file__).resolve().parents[1]
dest = root / 'examples/cds'
source = dest / 'NC_000932.gb'
expected_sha256 = 'a8b5d8239001f56a5b8b3ff047b10338b839329cf594aad36bfa4755a0dfb480'
if hashlib.sha256(source.read_bytes()).hexdigest() != expected_sha256:
    raise SystemExit('Original GenBank source hash changed; investigate before rebuilding')
record = SeqIO.read(source, 'genbank')
selected = []
for feature in record.features:
    if feature.type != 'CDS' or feature.qualifiers.get('gene', [''])[0] not in ('rps16', 'psbK', 'petD'):
        continue
    selected.append(dict(feature_id=feature.qualifiers['protein_id'][0],
        gene=feature.qualifiers['gene'][0], accver=record.id,
        strand='+' if feature.location.strand == 1 else '-',
        codon_start=int(feature.qualifiers['codon_start'][0]),
        genetic_code=int(feature.qualifiers['transl_table'][0]),
        segments=[dict(start=int(p.start)+1, end=int(p.end), segment_order=i+1)
                  for i, p in enumerate(feature.location.parts)],
        expected_protein=feature.qualifiers['translation'][0],
        source_location=str(feature.location), qualifiers=feature.qualifiers))
text = str(record.seq)
fasta = '>'+record.id+' '+record.description+'\n'+'\n'.join(text[i:i+70] for i in range(0,len(text),70))+'\n'
(dest/'NC_000932.fna').write_text(fasta)
(dest/'expected.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    fasta_sha256=hashlib.sha256(fasta.encode()).hexdigest(), parser='Biopython '+Bio.__version__,
    models=selected), indent=2)+'\n')
print('Extracted',len(selected),'annotated CDS and unchanged forward sequence')
