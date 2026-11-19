"""Save a derived handbook table and its translations from reviewed bilingual inputs."""
import json
import runpy
from . import ROOT


def write_page(slug, bodies, check=False):
    web = ROOT / 'website'
    sync = runpy.run_path(str(ROOT/'tools/sync-website-languages.py'))
    source = bodies['zh']
    canonical = web/'zh'/(slug+'.md')
    saved = web/'translations/en'/(slug+'.json')
    if check:
        if canonical.read_text() != source:
            raise ValueError('Stale handbook table: ' + slug)
        # Check the saved translation against the same current source contracts.
        translated = sync['translated_page'](source, json.loads(saved.read_text()))
        if translated != bodies['en']:
            raise ValueError('Stale English handbook table: ' + slug)
        return
    entries = {}
    def paragraphs(text):
        return [(kind, block) for kind, block in sync['blocks'](text) if kind != 'space']
    han, english = paragraphs(source), paragraphs(bodies['en'])
    if len(han) != len(english):
        raise ValueError('Unmatched table blocks: ' + slug)
    for index, ((kind, block), (en_kind, en)) in enumerate(zip(han, english)):
        if kind != en_kind:
            raise ValueError(f'{slug}/{index}: translated structure differs')
        if kind in ('code', 'directive', 'html'):
            continue
        entries[str(index)] = dict(source_sha256=sync['sha'](block), status='reviewed',
                                   translation=sync['mask_target'](block, en))
    record = dict(source_sha256=sync['sha'](source), blocks=entries)
    sync['translated_page'](source, record)
    canonical.write_text(source)
    saved.write_text(json.dumps(record, ensure_ascii=False, indent=2)+'\n')
