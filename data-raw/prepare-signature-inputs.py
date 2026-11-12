#!/usr/bin/env python3
"""Re-encode unchanged cached annotations, never run a signature analysis.

Missing MD5/status/run-date/description fields stay '-'. These TSVs are explicit
format round trips from cached annotations, not claimed InterProScan outputs.
"""
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
folder=root/'inst/extdata'
manifest_path=folder/'signature-sources.json'
items=json.loads(manifest_path.read_text())
for item in items:
    path=folder/item['file']
    if hashlib.sha256(path.read_bytes()).hexdigest()!=item['sha256']:
        raise SystemExit('Cached annotation bytes changed: '+item['file'])
    data=json.loads(path.read_text())
    rows=[]
    for entry in data['results']:
        metadata=entry['metadata']
        for protein in entry['proteins']:
            for location in protein['entry_protein_locations']:
                for fragment in location['fragments']:
                    score=location.get('score')
                    fields=[protein['accession'],'-',str(protein['protein_length']),
                        metadata['source_database'],metadata['accession'],metadata['name'],
                        str(fragment['start']),str(fragment['end']),
                        '-' if score is None else str(score),'-','-',
                        metadata.get('integrated') or '-','-']
                    if any(any(c in field for c in '\t\r\n') for field in fields):
                        raise SystemExit('Metadata cannot be re-encoded without editing a TSV field')
                    rows.append('\t'.join(fields))
    output=folder/path.name.replace('-pfam.json','-signatures.tsv')
    header=['# Re-encoding of the companion unchanged cached entry JSON.',
        '# Existing AA annotations only; no scanner or domain analysis was run.',
        '# Missing MD5, status, run date and integrated description remain -.']
    output.write_text('\n'.join(header+rows)+'\n')
    item['derived_tsv']={'file':output.name,'sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
        'rule':'One row per original location fragment; preserve accession case, AA boundaries, library/signature and reported score. Missing metadata remains -. This is not a new scan output.',
        'source_matches':sum(len(p['entry_protein_locations']) for e in data['results'] for p in e['proteins']),
        'output_fragments':len(rows)}
    print(output.name,len(rows),'original fragments re-encoded; no analysis')
manifest_path.write_text(json.dumps(items,indent=2)+'\n')
