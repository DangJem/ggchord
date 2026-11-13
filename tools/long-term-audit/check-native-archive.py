"""Verify actual local native artifacts; no rendering, downloads or release checks."""
from pathlib import Path
from zipfile import ZipFile
import glob
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[2]
digest = lambda value: hashlib.sha256(value).hexdigest()
snapshot_path = ROOT/'.site/source-snapshot.json'
snapshot = json.loads(snapshot_path.read_text())
assert snapshot['status'] == 'verified unchanged across rendering and finalization'
for name, expected in snapshot['render_inputs'].items():
    assert digest((ROOT/name).read_bytes()) == expected, name
engine = {str(p.relative_to(ROOT)): digest(p.read_bytes())
          for p in sorted((ROOT/'R').glob('*.R'))}
assert engine == {k: v for k, v in snapshot['render_inputs'].items()
                  if k.startswith('R/') and k.endswith('.R')}

scripts = sorted((ROOT/'.site/assets/examples').glob('*-en.R'))
assert scripts
literals = sorted({item for p in scripts for item in
                   re.findall(r'''["'](examples/[^"'\n]+)["']''', p.read_text())})
resolved = set()
for item in literals:
    matches = [Path(p) for p in glob.glob(str(ROOT/item)) if Path(p).is_file()]
    assert matches, item
    resolved.update(str(p.relative_to(ROOT)) for p in matches)
with ZipFile(ROOT/'.site/assets/tutorial-inputs.zip') as archive:
    hashes = json.loads(archive.read('input-hashes.json'))
    assert hashes == snapshot['archive_inputs']
    for name, expected in hashes.items():
        assert digest(archive.read(name)) == digest((ROOT/name).read_bytes()) == expected, name
    for name in resolved:
        assert archive.read(name) == (ROOT/name).read_bytes(), name
    shared_source = (ROOT/'tools/documentation-examples.R').read_bytes()
    assert archive.read('tools/documentation-examples.R') == shared_source
    assert (ROOT/'.site/assets/documentation-examples.R').read_bytes() == shared_source

result = dict(status='passed', scope='Actual native source snapshot, ZIP inputs and expanded English R references',
              archive_inputs_verified=len(hashes), expanded_english_scripts=len(scripts),
              literal_code_paths=len(literals), resolved_example_code_files=len(resolved),
              source_snapshot_sha256=digest(snapshot_path.read_bytes()), runtime_sha256=engine,
              verifier_sha256=digest(Path(__file__).read_bytes()),
              not_checked=['native full link/search checker', 'browser and plot quality',
                           'release package/device/platform checks'])
(ROOT/'tools/long-term-audit/v0.27.0-native-archive-results.json').write_text(
    json.dumps(result, indent=2)+'\n')
print(f'PASS: {len(hashes)} archive inputs and {len(resolved)} expanded code input files; source snapshot verified')
