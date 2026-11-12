#!/usr/bin/env python3
"""Verify actual downloadable archive bytes, not just a source-file inventory."""
import hashlib
import json
import re
import zipfile
from pathlib import Path

root = Path(__file__).resolve().parents[2]
site = root / '.site'
config = (root / '_pkgdown.yml').read_text()
topics = [re.search(r'\\name\{([^}]+)\}', p.read_text())[1] for p in (root / 'man').glob('*.Rd')]
assert all('`' + topic + '`' in config for topic in topics)
archive_path = site / 'assets/tutorial-inputs.zip'
with zipfile.ZipFile(archive_path) as archive:
    hashes = json.loads(archive.read('input-hashes.json'))
    assert hashes == json.loads((site / 'assets/input-hashes.json').read_text())
    for path, sha in hashes.items():
        assert hashlib.sha256(archive.read(path)).hexdigest() == sha
        assert hashlib.sha256((root / path).read_bytes()).hexdigest() == sha
    references = set()
    for script in (site / 'assets/examples').glob('*.R'):
        for pattern in re.findall(r'''["'](examples/[^"'\n]+)["']''', script.read_text()):
            matches = list(root.glob(pattern))
            assert matches, pattern
            references.update(str(path.relative_to(root)) for path in matches if path.is_file())
    assert references <= set(hashes), references - set(hashes)
    historical_present = 'examples/genome-alignments/nucmer.delta' in archive.namelist()
    assert not historical_present
result = dict(scope='Actual generated archive bytes and complete downloadable script file references; no deployment or complete site claim',
    inputs=len(hashes), literal_file_references=len(references), all_bytes_match_source=True,
    all_script_inputs_present=True, historical_nucleotide_data_in_archive=historical_present,
    native_topics_covered=len(topics), archive_sha256=hashlib.sha256(archive_path.read_bytes()).hexdigest())
(Path(__file__).parent / 'v0.27.0-site-input-results.json').write_text(json.dumps(result, indent=2) + '\n')
print('Verified archive:', len(hashes), 'source hashes and', len(references), 'script file references')
