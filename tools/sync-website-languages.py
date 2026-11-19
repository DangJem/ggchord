#!/usr/bin/env python3
"""Generate Han Markdown and bilingual metadata from English and reviewed translations."""
import argparse
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WEB = ROOT / 'website'
CONFIG = WEB / 'content-policy.json'
FENCE = re.compile(r'^\s*(`{3,}|~{3,})')
LIST = re.compile(r'^\s*(?:[-*+] |\d+[.)] )')
SHARED = {'space', 'code', 'directive', 'html'}
TOKEN = re.compile(r'`[^`\n]+`|\{#[^}\n]+\}|\{\{[^}\n]+\}\}|<[^>\n]+>|(?<=\]\()[^)\n]*(?=\))')
# The atlas formats its missing-default explanations as inline code.
ATLAS_LABELS = {
    '`No formal default; requiredness follows the contract`': '`无形式默认；是否必需见契约`',
    '`No formal default; see upstream contract`': '`无形式默认；见上游契约`',
}


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def json_text(value):
    return json.dumps(value, ensure_ascii=False, indent=2) + '\n'


def digest_json(value):
    return sha(json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':')))


def write_if_changed(path, text):
    if not path.exists() or path.read_text() != text:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)


def block_kind(line):
    if not line.strip():
        return 'space'
    if FENCE.match(line):
        return 'code'
    if re.fullmatch(r'\s*\{\{[^\n]+\}\}\s*', line):
        return 'directive'
    if re.match(r'^#{1,6} ', line):
        return 'heading'
    if line.lstrip().startswith('|'):
        return 'table'
    if LIST.match(line):
        return 'list'
    if re.fullmatch(r'(?:\s*<[^>]+>)+\s*', line):
        return 'html'
    return 'paragraph'


def blocks(text):
    lines = text.splitlines(keepends=True)
    start = 0
    while start < len(lines):
        end = start + 1
        kind = block_kind(lines[start])
        if kind == 'code':
            fence = FENCE.match(lines[start])[1]
            close = re.compile(r'^\s*' + fence[0] + '{' + str(len(fence)) + r',}\s*$')
            while end < len(lines) and not close.match(lines[end]):
                end += 1
            if end == len(lines):
                raise ValueError(f'Unclosed code fence at line {start + 1}')
            end += 1
        elif kind == 'space':
            while end < len(lines) and not lines[end].strip():
                end += 1
        elif kind in ('paragraph', 'list'):
            while end < len(lines) and block_kind(lines[end]) == 'paragraph' and not lines[end].lstrip().startswith('<'):
                end += 1
        yield kind, ''.join(lines[start:end])
        start = end


def localized_token(token):
    if token in ATLAS_LABELS:
        return ATLAS_LABELS[token]
    return re.sub(r'(^|/)MIGRATION\.md(?=$|#)', r'\1MIGRATION-Hans.md', token)


def validate_tokens(source, target):
    expected = Counter(localized_token(token) for token in TOKEN.findall(source))
    actual = Counter(TOKEN.findall(target))
    # Existing prose occasionally prints an identifier without its backticks.
    # Count literal, delimited occurrences outside all code spans, links and HTML.
    bare_prose = TOKEN.sub(' ', target)
    missing = expected - actual
    for token, count in list(missing.items()):
        if token.startswith('`'):
            name = re.escape(token[1:-1])
            bare = re.findall(r'(?<![a-zA-Z0-9_.])' + name + r'(?![a-zA-Z0-9_.])', bare_prose)
            if len(bare) >= count:
                del missing[token]
    if missing or actual - expected:
        raise ValueError(f'Shared code, links, markup or anchors differ: missing {list(missing)}, extra {list(actual - expected)}')


def validate_translation(source, target, kind=None):
    if not isinstance(target, str) or not target.strip():
        raise ValueError('Empty translation')
    if kind:
        if [k for k, _ in blocks(target)] != [kind]:
            raise ValueError(f'Translation changes the {kind} block structure')
        if kind == 'heading' and source.split(' ', 1)[0] != target.split(' ', 1)[0]:
            raise ValueError('Heading level differs')
        if kind == 'table' and len(re.findall(r'(?<!\\)\|', source)) != len(re.findall(r'(?<!\\)\|', target)):
            raise ValueError('Table column count differs')
    validate_tokens(source, target)


def english_blocks(source):
    return ((kind, text) for kind, text in blocks(source) if kind != 'space')


def checked_entry(entries, key, text, kind=None):
    entry = entries.get(key, {})
    if entry.get('source_sha256') != sha(text) or entry.get('status') != 'reviewed':
        raise ValueError(f'Missing, stale or unreviewed translation at {key}')
    try:
        translation = entry['translation']
        validate_translation(text, translation, kind)
    except (ValueError, KeyError) as error:
        raise ValueError(f'{key}: {error}') from error
    return translation


def reviewed_entry(source, target):
    return dict(source_sha256=sha(source), status='reviewed', translation=target)


def reviewed_page(source, target):
    originals = list(english_blocks(source))
    translations = list(english_blocks(target))
    if len(originals) != len(translations):
        raise ValueError(f'Block counts differ: {len(originals)} English, {len(translations)} Han')
    entries = {}
    for index, ((kind, text), (translated_kind, translation)) in enumerate(zip(originals, translations)):
        if translated_kind != kind:
            raise ValueError(f'Block {index}: expected {kind}, found {translated_kind}')
        if kind in SHARED:
            if text != translation:
                raise ValueError(f'Block {index}: shared {kind} differs')
        else:
            try:
                validate_translation(text, translation, kind)
            except ValueError as error:
                raise ValueError(f'Block {index}: {error}') from error
            entries[str(index)] = reviewed_entry(text, translation)
    return dict(source_sha256=sha(source), blocks=entries)


def translated_page(source, saved):
    if saved.get('source_sha256') != sha(source):
        raise ValueError('English source changed; review the Han translation')
    entries = saved['blocks']
    output = []
    used = set()
    index = 0
    for kind, text in blocks(source):
        if kind in SHARED:
            output.append(text)
        else:
            key = str(index)
            output.append(checked_entry(entries, key, text, kind))
            used.add(key)
        if kind != 'space':
            index += 1
    if used != set(entries):
        raise ValueError('Orphan translation blocks')
    return ''.join(output)


def map_strings(value, fn, path=()):
    if isinstance(value, dict):
        if set(value) == {'en'} and isinstance(value['en'], str):
            return fn('/'.join(map(str, path)), value['en'])
        if 'zh' in value and isinstance(value['zh'], str):
            raise ValueError(f'Canonical structured content contains a Han string at {path}')
        return {k: map_strings(v, fn, (*path, k)) for k, v in value.items()}
    if isinstance(value, list):
        return [map_strings(v, fn, (*path, i)) for i, v in enumerate(value)]
    return value


def split_content(paired):
    """Extract a canonical English tree and its explicitly supplied Han strings."""
    strings = {}

    def visit(value, path):
        if isinstance(value, dict):
            if set(value) == {'en', 'zh'}:
                text, target = value['en'], value['zh']
                key = '/'.join(map(str, path))
                validate_translation(text, target)
                strings[key] = reviewed_entry(text, target)
                return {'en': text}
            if 'en' in value or 'zh' in value:
                raise ValueError(f'Incomplete language pair at {path}')
            return {k: visit(v, (*path, k)) for k, v in value.items()}
        if isinstance(value, list):
            return [visit(v, (*path, i)) for i, v in enumerate(value)]
        return value

    source = visit(paired, ())
    return source, dict(source_sha256=digest_json(source), strings=strings)


def resolve_content(source, saved):
    if saved.get('source_sha256') != digest_json(source):
        raise ValueError('Structured English source changed; review the Han translation')
    entries = saved['strings']
    used = set()

    def pair(key, text):
        used.add(key)
        return {'zh': checked_entry(entries, key, text), 'en': text}

    result = map_strings(source, pair)
    if used != set(entries):
        raise ValueError('Orphan structured translations')
    return result


def generate(row, kind):
    source = (WEB / row['source']).read_text()
    saved = json.loads((WEB / row['translations']).read_text())
    if kind == 'markdown':
        return translated_page(source, saved), len(saved['blocks'])
    return json_text(resolve_content(json.loads(source), saved)), len(saved['strings'])


def policy():
    settings = json.loads(CONFIG.read_text())
    if (settings['source_language'], settings['generated_language']) != ('en', 'zh'):
        raise ValueError('The website content policy requires an English source and a Han translation')
    return settings


def check_inventory(settings):
    markdown = settings['markdown']
    structured = settings['structured']
    all_rows = markdown + structured
    for key in ('source', 'translations', 'output'):
        paths = [r[key] for r in all_rows]
        if len(paths) != len(set(paths)):
            raise ValueError(f'Duplicate {key} path in website/content-policy.json')
    for pattern, expected in [
        ('en/*.md', {r['source'] for r in markdown}),
        ('zh/*.md', {r['output'] for r in markdown}),
        ('content-en/**/*.json', {r['source'] for r in structured}),
        ('translations/zh/**/*.json', {r['translations'] for r in all_rows}),
    ]:
        actual = {p.relative_to(WEB).as_posix() for p in WEB.glob(pattern)}
        unregistered = sorted(actual - expected)
        if unregistered:
            raise ValueError(f'Files missing from website/content-policy.json: {", ".join(unregistered)}')


def sync(check=True):
    settings = policy()
    check_inventory(settings)
    counts = dict(pages=0, extra_markdown=0, structured_files=0, translated_blocks=0, translated_strings=0)
    outputs = {}
    for kind in ('markdown', 'structured'):
        for row in settings[kind]:
            try:
                output, count = generate(row, kind)
            except (OSError, ValueError, KeyError) as error:
                raise ValueError(f"{row['source']} + {row['translations']}: {error}") from error
            path = WEB / row['output']
            if check and (not path.exists() or path.read_text() != output):
                raise ValueError(f"Generated output differs: {row['output']}. Run tools/sync-website-languages.py")
            outputs[path] = output
            if kind == 'markdown':
                counts['pages' if row['source'].startswith('en/') else 'extra_markdown'] += 1
                counts['translated_blocks'] += count
            else:
                counts['structured_files'] += 1
                counts['translated_strings'] += count
    if not check:
        for path, text in outputs.items():
            write_if_changed(path, text)
    return counts


def prepare():
    """Only an unchanged English source string may keep its existing review."""
    settings = policy()
    drafts = 0
    for kind, field in [('markdown', 'blocks'), ('structured', 'strings')]:
        for row in settings[kind]:
            source = (WEB / row['source']).read_text()
            path = WEB / row['translations']
            previous = json.loads(path.read_text())[field] if path.exists() else {}
            by_hash = {entry['source_sha256']: entry for entry in previous.values()}
            entries = {}

            def keep_or_draft(key, text):
                nonlocal drafts
                digest = sha(text)
                old = previous.get(key, {})
                if old.get('source_sha256') == digest:
                    entry = old
                elif digest in by_hash and kind == 'markdown':
                    entry = by_hash[digest]
                else:
                    entry = dict(source_sha256=digest, status='needs_review', source_for_review=text,
                                 translation=old.get('translation', ''))
                    drafts += 1
                entries[key] = entry

            if kind == 'markdown':
                digest = sha(source)
                for index, (block_type, text) in enumerate(english_blocks(source)):
                    if block_type not in SHARED:
                        keep_or_draft(str(index), text)
            else:
                content = json.loads(source)
                digest = digest_json(content)
                map_strings(content, keep_or_draft)
            write_if_changed(path, json_text({'source_sha256': digest, field: entries}))
    return drafts


def review_outputs(row, kind, draft_path):
    output_path = WEB / row['output']
    if Path(draft_path).resolve() == output_path.resolve():
        raise ValueError('Supply a separately reviewed translation draft, not the generated file')
    source = (WEB / row['source']).read_text()
    draft = Path(draft_path).read_text()
    if kind == 'markdown':
        record = reviewed_page(source, draft)
        output = translated_page(source, record)
    else:
        expected, record = split_content(json.loads(draft))
        if expected != json.loads(source):
            raise ValueError('Reviewed paired JSON differs from the English source or its shared metadata')
        output = json_text(resolve_content(expected, record))
    return {WEB / row['translations']: json_text(record), output_path: output}


def record_review(kind, name, draft_path):
    field = 'source' if kind == 'markdown' else 'output'
    row = next((r for r in policy()[kind] if r[field] == name), None)
    if not row:
        raise ValueError(f'{name}: register the source and output in website/content-policy.json first')
    for path, content in review_outputs(row, kind, draft_path).items():
        write_if_changed(path, content)


def review_directory(folder):
    settings = policy()
    rows = {r['source']: r for r in settings['markdown']}
    drafts = sorted(Path(folder).glob('*.md'))
    if not drafts:
        raise ValueError(f'No reviewed Markdown drafts in {folder}')
    outputs = {}
    for draft in drafts:
        source = 'en/' + draft.name
        if source not in rows:
            row = dict(source=source, translations=f'translations/zh/{draft.stem}.json',
                       output='zh/' + draft.name)
            settings['markdown'].append(row)
        else:
            row = rows[source]
        try:
            outputs.update(review_outputs(row, 'markdown', draft))
        except (OSError, ValueError) as error:
            raise ValueError(f'{draft}: {error}') from error
    # Validate every supplied pair before writing anything, including new registrations.
    for path, content in outputs.items():
        write_if_changed(path, content)
    write_if_changed(CONFIG, json_text(settings))
    return len(drafts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', action='store_true')
    group.add_argument('--prepare', action='store_true', help='Prepare changed English prose for translation')
    group.add_argument('--reviewed-page', nargs=2, metavar=('ENGLISH_SOURCE', 'REVIEWED_HAN_MARKDOWN'))
    group.add_argument('--reviewed-content', nargs=2, metavar=('OUTPUT_NAME', 'REVIEWED_PAIRED_JSON'))
    group.add_argument('--reviewed-directory', metavar='REVIEWED_HAN_DIRECTORY',
                       help='Register reviewed Han drafts matching website/en/*.md and write their generated outputs')
    args = parser.parse_args()
    try:
        if args.prepare:
            print('Prepared', prepare(), 'new or changed translation drafts')
        elif args.reviewed_directory:
            print('Saved reviewed Han translations for', review_directory(args.reviewed_directory), 'pages')
        elif args.reviewed_page or args.reviewed_content:
            kind = 'markdown' if args.reviewed_page else 'structured'
            name, path = args.reviewed_page or args.reviewed_content
            record_review(kind, name, path)
            print('Saved reviewed Han translation for', name)
        else:
            counts = sync(check=args.check)
            print('PASS: English source → reviewed Han; ' + ', '.join(f'{v} {k}' for k, v in counts.items()))
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f'{error}\n')


if __name__ == '__main__':
    main()
