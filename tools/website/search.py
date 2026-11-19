"""Structured local search index and downloadable examples."""
import html,json,re,hashlib,shutil
from urllib.parse import quote
from website import ROOT
from website.handbook import EXAMPLE, WEB, manifest
from website.registry import page_kind, pages

def complete_examples(assets):
    folder = assets/'examples'
    folder.mkdir(exist_ok=True)
    specs = manifest()
    records = []
    for lang in ('en', 'zh'):
        for page in pages():
            slug = page['slug']
            if page.get('kind') == 'alias' or slug == 'reference':
                continue
            text = (WEB/lang/(slug + '.md')).read_text()
            for index, match in enumerate(EXAMPLE.finditer(text), 1):
                example_id = match[1]
                if example_id not in specs:
                    raise ValueError(f'{lang}/{slug}: unknown example {example_id}')
                source = WEB/'examples'/(example_id + '.R')
                name = f'{slug}-{index}-{lang}.R'
                shutil.copyfile(source, folder/name)
                records.append(dict(page=slug, index=index, language=lang, id=example_id,
                                    file='assets/examples/'+name,
                                    sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
    (assets/'examples.json').write_text(json.dumps(records, ensure_ascii=False, indent=2)+'\n')


def search_records(out,lang,page_records):
    records=[]
    for page in page_records:
        kind=page_kind(page['url'])
        if page['url']=='reference/index.html':
            page['kind']='tool'
        elif page['url']=='faq.html':
            page['kind']='faq'
        elif kind=='api':
            page['kind']='function'
        elif kind in ('tutorial','guide','case','index','gallery'):
            page['kind']='guide'
        else:
            page['kind']='tool'
        page['keywords']=page['url'].removesuffix('.html')
        page['summary']=page['text'][:240]
        records.append(page)
        # Hit-level context uses only public rendered prose; executable blocks are excluded.
        path=out/lang/page['url']
        raw=path.read_text()
        main=re.search(r'<main\b[^>]*>(.*?)</main>',raw,re.S)
        body=main[1] if main else ''
        headings=list(re.finditer(r'<h([23])\b[^>]*id="([^"]+)"[^>]*>(.*?)</h\1>',body,re.S))
        if page['url']!='parameter-atlas.html':
            for i,h in enumerate(headings):
                text=body[h.end():headings[i+1].start() if i+1<len(headings) else len(body)]
                text=re.sub(r'<pre\b[^>]*>.*?</pre>','',text,flags=re.S)
                text=html.unescape(re.sub(r'<[^>]+>',' ',text))
                text=re.sub(r'\s+',' ',text).strip()
                section=html.unescape(re.sub(r'<[^>]+>','',h[3]))
                records.append(dict(kind=page['kind'],url=page['url']+'#'+h[2],title=page['title']+' · '+section,
                    section=section,keywords=section,text=text[:6000],summary=text[:300]))

    # Search the shipped signatures and the localized table cells readers see.
    # Grouped Rd topics can document several exported aliases.
    def plain(markup):
        return ' '.join(html.unescape(re.sub(r'<[^>]*>', ' ', markup)).split())

    exports = set(re.findall(r'^export\(([^)]+)\)', (ROOT/'NAMESPACE').read_text(), re.M))
    titles = {p['slug']: p['title'][lang] for p in pages()}
    for path in sorted((ROOT/'man').glob('*.Rd')):
        names = set(re.findall(r'\\alias\{([^}]+)\}', path.read_text())) & exports
        if not names:
            continue
        url = 'reference/' + path.stem + '.html'
        raw = (out/lang/url).read_text()
        title = plain(re.search(r'<h1\b[^>]*>(.*?)</h1>', raw, re.S)[1])
        description = re.search(r'<div class="ref-description[^"]*">(.*?)</div>', raw, re.S)
        summary = plain(description[1]) if description else title
        target = url + ('#ref-usage' if 'id="ref-usage"' in raw else '')
        for api in sorted(names):
            records.append(dict(kind='function', url=target, title=api, section=title,
                                keywords=api + ' ' + title, text=summary, summary=summary))
    atlas = (out/lang/'parameter-atlas.html').read_text()
    for row in re.finditer(r'<tr\b[^>]* id="(param-[^"]+)"[^>]*>(.*?)</tr>', atlas, re.S):
        anchor = html.unescape(row[1])
        api, parameter = anchor.removeprefix('param-').split('-', 1)
        cells = [plain(cell) for cell in re.findall(r'<td\b[^>]*>(.*?)</td>', row[2], re.S)]
        if len(cells) != 3:
            raise ValueError(f'{lang}/parameter-atlas: expected parameter, function and usage in {anchor}')
        summary = ' · '.join(cell for cell in cells if cell)
        records.append(dict(kind='parameter', url='parameter-atlas.html#' + quote(anchor, safe='._-'),
                            title=api + ' · ' + parameter, section=titles['parameter-atlas'],
                            keywords=api + ' ' + parameter, text=' '.join(cells), summary=summary))
    for intent in json.loads((WEB/'search-intents.json').read_text()):
        records.append(dict(kind='task',url=intent['url'],title=intent['title'][lang],queries=intent['queries'],keywords=' '.join(intent['queries']),text=intent['title'][lang],summary=intent['summary'][lang]))
    return records
