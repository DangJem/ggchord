#!/usr/bin/env python3
"""Retry failed native jobs from unchanged inputs, retaining every initial result."""
import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
spec = importlib.util.spec_from_file_location('native_builder', ROOT/'tools/build-website.py')
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--record', default='tools/site-audit/manual-build-recovery.json',
                    help='Repository-relative JSON result path under tools/site-audit')
parser.add_argument('--language', choices=('en', 'zh'),
                    help='Retry one failed locale while the other locale continues')
parser.add_argument('--defer-finalize', action='store_true',
                    help='Retain retry results without finalizing an active build')
args = parser.parse_args()
record = (ROOT/args.record).resolve()
assert record.parent == (ROOT/'tools/site-audit').resolve() and record.suffix == '.json'
out = ROOT/'.site'
snapshot = json.loads((out/'source-snapshot.json').read_text())
assert snapshot['status'] == 'started; completion unverified'
process_file = out/'render-processes.json'
if process_file.exists():
    processes = json.loads(process_file.read_text())
else:
    assert args.language and args.defer_finalize, 'Main build is still active'
    workers = json.loads((out/(args.language+'-workers.json')).read_text())
    assert workers[-1]['return_code'] != 0
    processes = dict(status='other locale still active', return_codes={args.language: 1})
assert any(processes['return_codes'].values()), 'No failed native build to recover'
payload = json.loads((out/'page-content.json').read_text())
builder.assert_render_source_snapshot(snapshot['render_inputs'], snapshot['archive_inputs'])
recovery = dict(initial_processes=processes, retries=[],
                scope='Native HTML retries only; exact original render/archive inputs required; no browser, package or publication checks')
if record.exists():
    previous = json.loads(record.read_text())
    assert previous['render_inputs'] == snapshot['render_inputs'], 'Retry source snapshot changed'
    recovery['retries'] = previous['retries']
recovery['render_inputs'] = snapshot['render_inputs']
for lang in ('en', 'zh'):
    if args.language and lang != args.language:
        continue
    path = out/(lang+'-workers.json')
    workers = json.loads(path.read_text())
    if processes['return_codes'][lang] == 0:
        continue
    # The original batch runner stops at the first failed job. Retry that batch
    # one article per process, then any not-yet-attempted remaining articles.
    failed = workers[-1]
    assert failed['return_code'] != 0
    if '--part=reference' in failed['options']:
        reference_start = int(next((x.split('=')[1] for x in failed['options']
                                   if x.startswith('--offset=')), '0'))
        jobs = [['--part=reference', '--offset='+str(i), '--count=1']
                for i in range(reference_start, len(list((ROOT/'man').glob('*.Rd'))))]
        start = 0
    else:
        jobs = []
        start = int(next(x.split('=')[1] for x in failed['options'] if x.startswith('--offset=')))
    jobs += [['--part=articles', '--offset='+str(i), '--count=1']
             for i in range(start, len(payload[lang]))]
    for index, options in enumerate(jobs):
        builder.assert_render_source_snapshot(snapshot['render_inputs'], snapshot['archive_inputs'])
        completed = [r for r in recovery['retries'] if r['language'] == lang and r['options'] == options
                     and r['return_code'] == 0]
        if completed:
            assert (out/completed[-1]['log']).exists()
            continue
        log = out/(lang+'-recovery-'+str(index)+'.log')
        with log.open('w') as stream:
            run = subprocess.run(['Rscript', '--vanilla', 'tools/render-website.R', lang, *options],
                                 cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
        attempt = dict(language=lang, options=options, return_code=run.returncode, log=log.name)
        recovery['retries'].append(attempt)
        record.write_text(json.dumps(recovery, indent=2)+'\n')
        if run.returncode:
            raise SystemExit('Retry failed: '+log.name+'\n'+log.read_text()[-5000:])
        builder.assert_render_source_snapshot(snapshot['render_inputs'], snapshot['archive_inputs'])
    failed['recovery'] = [r for r in recovery['retries'] if r['language'] == lang]
    path.write_text(json.dumps(workers, indent=2)+'\n')
if args.defer_finalize:
    recovery['status'] = 'retry complete; finalization deferred'
    record.write_text(json.dumps(recovery, indent=2)+'\n')
    print('PASS: selected failed locale recovered; original worker results retained')
    raise SystemExit(0)
processes['initial_return_codes'] = processes['return_codes']
processes['return_codes'] = dict(en=0, zh=0)
processes['recovery_record'] = str(record.relative_to(ROOT))
(out/'render-processes.json').write_text(json.dumps(processes, indent=2)+'\n')
builder.assert_render_source_snapshot(snapshot['render_inputs'], snapshot['archive_inputs'])
builder.audit_reference_output(out)
builder.finalize(out, argparse.Namespace(publication=False), snapshot['render_inputs'], snapshot['archive_inputs'])
recovery['status'] = 'recovered and finalized with unchanged source inputs'
recovery['initial_workers_preserved'] = True
record.write_text(json.dumps(recovery, indent=2)+'\n')
manifest_path = out/'build-manifest.json'
manifest = json.loads(manifest_path.read_text())
manifest['source_sha256'][str(record.relative_to(ROOT))] = hashlib.sha256(record.read_bytes()).hexdigest()
manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n')
print('PASS: failed native jobs recovered and finalized; initial failures retained')
