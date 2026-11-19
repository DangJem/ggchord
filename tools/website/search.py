"""Structured local search index and downloadable examples."""
import html,json,re,hashlib
from pathlib import Path
from urllib.parse import quote
from documentation_lib import ROOT,WEB,PAGES,RECIPES,snippets,expand
from website.registry import page_kind

def complete_examples(assets):
    folder=assets/'examples';folder.mkdir(exist_ok=True)
    records = []
    for lang in ('en', 'zh'):
        for page in PAGES:
            slug = page['slug']
            if page.get('kind') == 'alias':
                continue
            text = (WEB/lang/(slug+'.md')).read_text()
            if not re.search(r'\{\{(?:code|prepare):|^```r', text, re.M):
                continue
            expanded = expand(text, lang, figures='../assets')
            blocks = re.findall(r'^```\{[^}]*\brun-example\b[^}]*\}\n(.*?)\n```', expanded, re.S | re.M)
            for index, block in enumerate(blocks, 1):
                code = block + '\n'
                name = f'{slug}-{index}-{lang}.R'
                (folder/name).write_text(code)
                records.append(dict(page=slug, index=index, language=lang, file='assets/examples/'+name,
                                    sha256=hashlib.sha256(code.encode()).hexdigest()))
    (assets/'examples.json').write_text(json.dumps(records,ensure_ascii=False,indent=2)+'\n')

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
    titles = {p['slug']: p['title'][lang] for p in PAGES}
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
        default = cells[1]
        summary = ('默认值 ' if lang == 'zh' else 'Default ') + default + ' · ' + cells[-1]
        records.append(dict(kind='parameter', url='parameter-atlas.html#' + quote(anchor, safe='._-'),
                            title=api + ' · ' + parameter, section=titles['parameter-atlas'],
                            keywords=api + ' ' + parameter, text=' '.join(cells), summary=summary))
    for r in json.loads((WEB/'teaching-tasks.json').read_text()):
        slug, _, anchor = r['page'].partition('#')
        records.append(dict(kind='task', url=slug + '.html' + ('#' + anchor if anchor else ''),
                            title=r['task'][lang], section=titles[slug],
                            keywords=r['api'] + ' ' + r['parameters'], text=r['task'][lang],
                            summary=r['api'] + ' · ' + r['parameters'] + ' · ' + r['modes']))
    for intent in json.loads((WEB/'search-intents.json').read_text()):
        records.append(dict(kind='task',url=intent['url'],title=intent['title'][lang],queries=intent['queries'],keywords=' '.join(intent['queries']),text=intent['title'][lang],summary=intent['summary'][lang]))
    return records
