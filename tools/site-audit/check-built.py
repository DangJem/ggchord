#!/usr/bin/env python3
"""Audit the actual local handbook, downloads, images, native help and historical URLs."""
import argparse
import hashlib
import json
import re
import struct
import subprocess
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urljoin, urlsplit
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
WEB = ROOT/'website'
LANGUAGES = ('en','zh')
VOID = set('area base br col embed hr img input link meta param source track wbr'.split())
BASE = 'https://local.invalid/ggchord/'
MARKER = re.compile(r'\{\{example:([a-z0-9-]+)\}\}')
RETIRED = {'geom_gene','geom_feature_plasmid','geom_feature_label_repel','position_plasmid','position_track'}


class Document(HTMLParser):
    """A small HTML tree, preserving inline whitespace and syntax-highlighted text."""
    def __init__(self, path):
        super().__init__(convert_charrefs=True)
        self.root = ET.Element('document')
        self.stack = [self.root]
        self.feed(path.read_text())
        self.close()

    def handle_starttag(self, tag, attrs):
        element = ET.SubElement(self.stack[-1], tag, dict(attrs))
        if tag not in VOID:
            self.stack.append(element)

    def handle_startendtag(self, tag, attrs):
        ET.SubElement(self.stack[-1], tag, dict(attrs))

    def handle_endtag(self, tag):
        for index in range(len(self.stack)-1, 0, -1):
            if self.stack[index].tag == tag:
                del self.stack[index:]
                break

    def handle_data(self, content):
        elem = self.stack[-1]
        if len(elem):
            elem[-1].tail = (elem[-1].tail or '')+content
        else:
            elem.text = (elem.text or '')+content


def read(path):
    return json.loads(path.read_text())


def text(elem):
    return ''.join(elem.itertext()) if elem is not None else ''


def styled(elem, css, tag=None):
    return [n for n in elem.iter(tag) if css in n.get('class','').split()]


def target(page, url):
    """Resolve the same relative URL a browser would see on a project Pages site."""
    result = urlsplit(urljoin(BASE+page, url))
    if result.netloc != 'local.invalid' or result.scheme != 'https':
        return None
    path = unquote(result.path)
    if not path.startswith('/ggchord/'):
        raise ValueError(f'{page}: link outside the website: {url}')
    name = path[len('/ggchord/'):]
    return name+('index.html' if path.endswith('/') else ''), unquote(result.fragment)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clipboard(details):
    code = details.find('.//pre//code')
    if code is None:
        code = details.find('.//pre')
    if code is None:
        raise ValueError('A code disclosure has no preformatted content')
    lines = styled(code, 'r-in', 'span')
    return '\n'.join(map(text,lines)) if lines else text(code)


def scripts_for_comparison(rows):
    # R's parser discards localized comments and preserves string literals,
    # call arguments, signatures and the order of all executable expressions.
    code = '''
    x <- jsonlite::fromJSON(file("stdin"), simplifyVector=FALSE)
    parse_input <- function(s) parse(text=s, keep.source=FALSE)
    for (id in names(x)) {
      v <- x[[id]]
      if (!identical(parse_input(v[[1]]), parse_input(v[[2]]))) stop(id, ": EN and Han R code differ")
    }'''
    subprocess.run(['Rscript','--vanilla','-e',code], input=json.dumps(rows), text=True, check=True)


def run_redirects(scripts, cases):
    # Exercise the generated redirects themselves, including the historical hash.
    code = """
const {scripts,cases} = JSON.parse(require('fs').readFileSync(0,'utf8'));
const vm = require('vm'), base = 'https://local.invalid/ggchord/';
const compiled = new Map(Object.entries(scripts).map(([name,code])=>[name,new vm.Script(code)]));
for (const [start,end] of cases) {
  let url = new URL(start,base);
  for (let i=0;i<10;i++){
    const name=url.pathname.slice(new URL(base).pathname.length), script=compiled.get(name);
    if(!script) break;
    let next;
    script.runInNewContext({location:{pathname:url.pathname, hash:url.hash, replace:s=>{next=new URL(s,url)}},
        addEventListener:()=>{}, localStorage:{getItem:()=>null}});
    if(!next || next.href===url.href) break;
    url=next;
  }
  if(url.href!==new URL(end,base).href) throw Error(start+' => '+url.href+'; expected '+end);
}"""
    subprocess.run(['node','-e',code],input=json.dumps(dict(scripts=scripts,cases=cases)), text=True,check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--site', type=Path, default=ROOT/'.site')
    parser.add_argument('--report', type=Path, help='Write a JSON summary at this path')
    parser.add_argument('--reference-inputs', type=Path, help='Save the actual EN/Han reference clipboard inputs for independent replay')
    args = parser.parse_args()
    site = args.site.resolve()
    if not (site/'build-manifest.json').is_file():
        parser.error(f'{site} has no completed build manifest')
    manifest = read(site/'build-manifest.json')
    registry = read(WEB/'pages.json')
    pages = {p['slug']:p for p in registry if p['kind']!='alias'}
    alias_pages = {p['slug']+'.html' for p in registry if p['kind']=='alias'}
    navigation = read(WEB/'navigation.json')
    groups = [g for g in navigation if g.get('navbar')]
    specs = read(WEB/'handbook-examples.json')
    home = read(WEB/'home.json')['modules']
    previews = {row['id']:row for row in read(WEB/'home-previews.json')}
    full_cases = {row['id']:row['page'] for row in read(WEB/'manual-figures.json')}
    example_cases = {p['canonical_path']:MARKER.findall((WEB/'en'/(p['slug']+'.md')).read_text())
                     for p in pages.values() if p['kind']=='case' and p['slug'].startswith('case-')}
    flex = read(WEB/'flexible-parameters.json')
    parameters = read(site/'parameter-directory.json')
    formals = {(r['api'],r['parameter']):r for r in parameters}
    rd = {p.stem:re.findall(r'\\alias\{([^}]+)\}',p.read_text()) for p in (ROOT/'man').glob('*.Rd')}
    exports = set(re.findall(r'^export\(([^)]+)\)',(ROOT/'NAMESPACE').read_text(),re.M))
    records = read(site/'assets/examples.json')
    downloads = {(r['language'],r['page'],r['index']):r for r in records}
    rendered = read(site/'assets/handbook/rendered.json')
    output = read(site/'reference-output.json')
    output_rows = {lang:{r['topic']:r for r in output['topics'][lang]} for lang in LANGUAGES}
    migrations = read(WEB/'page-migrations.json')
    routes = read(site/'legacy-routes.json')
    retired = set(read(site/'retired-routes.json'))
    identifiers, scripts, references, signatures = {}, {}, {}, {}
    links, cases, errors = set(), set(), []
    used_downloads = set()
    counts = Counter()
    figures = {}
    body_anchors = {}
    expected_pages = {p['canonical_path'] for p in pages.values()}
    canonical_topics = {'reference/'+slug+'.html' for slug in rd}
    expected_pages.update(canonical_topics)
    expected_pages.update(('authors.html','LICENSE.html','LICENSE-text.html','404.html'))

    def require(ok, message):
        if not ok:
            errors.append(message)

    def destination(page, url):
        try:
            return target(page,url)
        except ValueError as error:
            errors.append(str(error))
            return None

    def remember(page, url):
        pair = destination(page,url)
        if pair:
            links.add((page,url,*pair))
        return pair

    require(len(formals)==len(parameters) and all(r['api'] in exports and r['api'] in rd.get(r['topic'],[]) for r in parameters),
            'Parameter directory names a duplicate or unknown exported argument')
    require(all((api,arg) in formals for api,args_here in flex.items() for arg in args_here),
            'Flexible argument missing from the runtime parameter directory')
    plots = {key for key,s in specs.items() if s['kind']=='plot'}
    require(set(rendered)==plots, 'Rendered plot inventory differs from registered plotting scripts')
    for name in sorted(plots):
        png = site/'assets/handbook'/(name+'.png')
        script = WEB/'examples'/(name+'.R')
        record = rendered.get(name,{})
        require(png.is_file() and record.get('png_sha256')==digest(png), name+': wrong rendered PNG digest')
        require(record.get('source_sha256')==digest(script), name+': plotted script differs from source')
        if png.is_file():
            with png.open('rb') as handle:
                header = handle.read(24)
            dims = tuple(round(specs[name][d]*specs[name].get('dpi',144)) for d in ('width','height'))
            require(len(header)==24 and header[:8]==b'\x89PNG\r\n\x1a\n' and
                    all(abs(x-y)<=1 for x,y in zip(struct.unpack('>II', header[16:24]), dims)),
                    name+': PNG format/dimensions differ from plotting manifest')
    require(output.get('status')=='passed', 'Native render reported failed example output')
    require(len(downloads)==len(records), 'Duplicate example download records')
    primary = ['index']+[g['landing'] for g in groups if g['landing']!='index']
    expected_nav = [pages[k]['canonical_path'] for k in primary]
    menus = {}
    for group in groups:
        keys = list(group['pages'])
        for child in navigation:
            if child.get('parent')==group['id']:
                keys.extend([child['landing'],*child['pages']])
        paths = [pages[key]['canonical_path'] for key in dict.fromkeys(keys) if key!=group['landing']]
        if paths:
            menus['site-menu-'+group['id']] = paths
    for path in sorted(site.rglob('*.html')):
        rel = path.relative_to(site).as_posix()
        parts = rel.split('/',1)
        lang, local = parts if len(parts)==2 and parts[0] in LANGUAGES else (None, rel)
        doc = Document(path).root
        all_ids = [n.get('id') for n in doc.iter() if n.get('id')]
        identifiers[rel] = set(all_ids)
        require(len(all_ids)==len(set(all_ids)),rel+': duplicate IDs')
        for node in doc.iter():
            for attr in ('href','src','data-code-url','data-search'):
                if node.get(attr):
                    remember(rel, node.get(attr))
            if node.get('srcset'):
                for item in node.get('srcset').split(','):
                    remember(rel,item.strip().split()[0])
            if node.tag=='img':
                require('alt' in node.attrib, rel+': image lacks an alt attribute')
        if lang:
            require(doc.find('.//html').get('lang')==('en' if lang=='en' else 'zh-Hans'),rel+': incorrect document language')
        main_content = doc.find('.//main')
        if rel in routes or rel in retired or (lang and local in migrations):
            inline = [text(n) for n in doc.iter('script') if 'location.replace' in text(n)]
            if inline:
                scripts[rel] = '\n'.join(inline)
        if lang and local in alias_pages:
            require(rel in scripts and not styled(doc,'code-disclosure') and not styled(doc,'site-primary-nav'),
                    rel+': historical article must be a redirect')
        if not lang or local not in expected_pages:
            continue
        require(main_content is not None,rel+': expected a canonical page, found a redirect or incomplete document')
        if main_content is None:
            continue
        main_ids = {n.get('id') for n in main_content.iter() if n.get('id')}
        body_anchors[rel] = main_ids
        require({'main','site-search','site-language'} <= identifiers[rel], rel+': missing accessibility, search or language entry')
        nav = styled(doc,'site-primary-nav')
        targets = [destination(rel,a.get('href','')) for a in styled(nav[0],'nav-link','a')] if nav else []
        require(targets==[(lang+'/'+n,'') for n in expected_nav],rel+': incorrect primary navigation')
        actual_menus = {n.get('id'):[destination(rel,a.get('href','')) for a in n.iter('a')]
                        for n in styled(doc,'site-nav-menu')}
        require(actual_menus=={key:[(lang+'/'+path,'') for path in values] for key,values in menus.items()},
                rel+': dropdown destinations differ from the navigation registry')
        other = 'zh' if lang=='en' else 'en'
        switch = doc.find('.//a[@id="site-language"]')
        require(switch is not None and destination(rel,switch.get('href',''))==(other+'/'+local,''),
                rel+': language switch does not lead to corresponding page')
        toggles = styled(doc,'site-nav-toggle','button')
        require({t.get('aria-controls') for t in toggles}==set(menus),rel+': missing menu control')
        for toggle in toggles:
            require(toggle.get('aria-controls') in identifiers[rel] and toggle.get('aria-expanded')=='false',
                    rel+': menu button has no matching controlled menu')
        for table in main_content.iter('table'):
            rows = table.findall('./tbody/tr')
            if not rows:
                continue
            for i in range(max(map(len,rows))):
                cells = [row[i] for row in rows if len(row)>i]
                values = [text(cell).strip() for cell in cells if text(cell).strip()]
                numeric = bool(values) and all(re.fullmatch(r'[+-]?(?:[0-9]+(?:[.,][0-9]+)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?%?',v) for v in values)
                headers = [r[i] for r in table.findall('./thead/tr') if len(r)>i]
                require(all(('numeric' in c.get('class','').split())==numeric for c in cells+headers),
                        f'{rel}: inconsistent table alignment in column {i+1}')
        if local in canonical_topics:
            topic = path.stem
            sections = {n.find('h2').get('id'):n for n in main_content.iter() if n.find('h2') is not None}
            example = sections.get('ref-examples')
            require(example is not None and 'ref-related' in main_ids,rel+': missing native examples or related entries')
            if example is None:
                continue
            blocks = styled(example,'code-disclosure','details')
            require(len(blocks)==1 and blocks[0].get('data-code-kind')=='r',rel+': expected one native R example')
            if not blocks:
                continue
            copied = clipboard(blocks[0])
            references[(lang,topic)] = copied
            if args.reference_inputs:
                folder = args.reference_inputs/lang
                folder.mkdir(parents=True, exist_ok=True)
                (folder/(topic+'.R')).write_text(copied+'\n')
            sigs = [clipboard(b) for b in styled(main_content,'code-disclosure','details') if b.get('data-code-kind')=='signature']
            signatures[(lang,topic)] = sigs
            if 'ref-usage' in main_ids:
                require(bool(sigs),rel+': copyable usage must be identified as a signature')
            outline = styled(doc,'page-nav')
            outline_targets = {a.get('href') for a in outline[0].iter('a')} if outline else set()
            require({'#'+key for key in sections if key} <= outline_targets,rel+': incomplete native page outline')
            imgs = [n.get('src') for n in example.iter('img')]
            figures[(lang,topic)] = imgs
            row = output_rows[lang].get(topic,{})
            require(imgs==row.get('figures') and (bool(imgs) or bool(styled(example,'r-out'))) and
                    not styled(example,'r-err'),rel+': missing native evaluated output or output audit mismatch')
            marked = set()
            for dl in main_content.iter('dl'):
                args_here = []
                for child in dl:
                    if child.tag=='dt' and child.get('id','').startswith('arg-'):
                        args_here = [v.strip() for v in text(child).split(',')]
                    if child.tag=='dd':
                        for badge in styled(child,'flexible-parameter'):
                            marked.update(args_here)
                            require(any(a.get('href','').startswith('../parameters.html#') for a in badge.iter('a')),
                                    rel+': flexible badge lacks a matching handbook section')
            require(marked==set(flex.get(topic,{})),f'{rel}: incorrect flexible parameter badges {marked}')
            if lang=='zh':
                require(bool(re.search('[\u3400-\u9fff]',text(main_content.find('.//h1')))) or
                        text(main_content.find('.//h1')).strip() in rd[topic], rel+': missing Han topic title')
            continue
        page = next((p for p in pages.values() if p['canonical_path']==local),None)
        if page is None:
            continue
        slug = page['slug']
        source = (WEB/lang/(slug+'.md')).read_text()
        expected = MARKER.findall(source)
        blocks = [n for n in main_content.iter() if n.get('data-example-kind')]
        dl_links = styled(main_content,'example-download','a')
        require(len(expected)==len(blocks)==len(dl_links),f'{rel}: marker/code/download counts {len(expected)}/{len(blocks)}/{len(dl_links)}')
        actual_figs = [destination(rel,img.get('src','')) for fig in styled(main_content,'handbook-output','figure') for img in fig.iter('img')]
        require(actual_figs==[('assets/handbook/'+name+'.png','') for name in expected if specs[name]['kind']=='plot'],
                rel+': missing, reordered or unrelated handbook plot output')
        for i,(name,block,anchor) in enumerate(zip(expected,blocks,dl_links),1):
            key = (lang,slug,i)
            used_downloads.add(key)
            registered = downloads.get(key,{})
            script = WEB/'examples'/(name+'.R')
            dest = destination(rel,anchor.get('href',''))
            require(registered.get('id')==name and dest==(registered.get('file'),''),rel+': incorrect download link for '+name)
            if dest is None:
                continue
            file = site/dest[0]
            copied = text(block.find('.//code'))
            require(block.get('data-example-kind')==specs[name]['kind'] and copied==script.read_text().rstrip('\n'),
                    rel+': displayed R code differs from '+name+'.R')
            require(file.is_file() and file.read_bytes()==script.read_bytes() and digest(file)==registered.get('sha256'),
                    rel+': downloaded bytes differ from '+name+'.R')
        if local=='index.html':
            cards = styled(main_content,'home-module','section')
            require(len(cards)==len(home)==4 and sum(len(g['figures']) for g in home)==7,rel+': expected four homepage groups and seven previews')
            for card,item in zip(cards,home):
                pictures = styled(card,'home-preview','figure')
                keys = [Path(img.get('src','')).stem for preview in pictures for img in preview.iter('img')]
                require(keys==item['figures'],rel+': wrong flagship preview in '+item['slug'])
                actions = [a.get('href') for a in card.iter('a') if {'module-entry','home-source'} & set(a.get('class','').split())]
                require(actions==[item['entry'],item['case']],rel+': try-it/complete case link mismatch')
                for key,preview in zip(keys,pictures):
                    links_to_figure = [a for a in preview.iter('a') if a.get('data-code-url')]
                    require(len(links_to_figure)==1,rel+': preview has no case and full figure links')
                    if not links_to_figure or key not in previews:
                        continue
                    case = full_cases[previews[key]['full']]
                    ids = example_cases.get(case,[])
                    anchor = links_to_figure[0]
                    require(anchor.get('data-code-url')==case and len(ids)==1 and
                            destination(rel,anchor.get('href',''))==('assets/handbook/'+ids[0]+'.png',''),
                            rel+': homepage expanded image does not match case script for '+key)
                    original = ROOT/previews[key]['file']
                    served = site/'assets'/original.name
                    require(served.is_file() and digest(served)==digest(original),rel+': homepage preview differs from source '+key)
        if '{{case-gallery}}' in source:
            linked_cases = []
            for card in styled(main_content,'gallery-card'):
                destinations = {a.get('href') for a in card.iter('a')}
                cases_here = destinations & set(example_cases)
                require(len(cases_here)==1 and len(destinations-cases_here)==1,rel+': gallery card lacks case/first-plot links')
                if len(cases_here)==1:
                    case = cases_here.pop()
                    linked_cases.append(case)
                    case_plots = example_cases[case]
                    images = [destination(rel,img.get('src','')) for img in card.iter('img')]
                    require(images==[('assets/handbook/'+key+'.png','') for key in case_plots],
                            rel+': gallery thumbnail differs from the complete example '+case)
            require(len(linked_cases)==len(set(linked_cases))==7 and set(linked_cases)==set(example_cases),
                    rel+': missing one of the seven public figure families')
        if '{{parameter-directory}}' in source:
            rows = {n.get('id'):n for n in main_content.iter('tr') if n.get('id','').startswith('param-')}
            require(set(rows)=={'param-'+api+'-'+arg for api,arg in formals},rel+': incomplete parameter directory')
            for api,arg in formals:
                row = rows.get('param-'+api+'-'+arg)
                if row is None:
                    continue
                entry = formals[(api,arg)]
                expected_link = 'reference/'+entry['topic']+'.html#'+entry['anchor']
                require(any(a.get('href')==expected_link for a in row.iter('a')),
                        f'{rel}: incorrect argument destination for {api}({arg})')
                supported = arg in flex.get(api,{})
                require(row.get('data-flex')==str(int(supported)) and bool(styled(row,'flexible-label','a'))==supported,
                        f'{rel}: incorrect flexible flag for {api}({arg})')
        if local=='reference/index.html':
            entries = {destination(rel,a.get('href','')) for group in styled(main_content,'reference-group')
                       for dl in group.iter('dl') for a in dl.iter('a')}
            require({(lang+'/'+key,'') for key in canonical_topics} <= entries,rel+': reference index omits a native help topic')
            require('reference-filter' in identifiers[rel],rel+': no reference filter')
        require(not re.search(r'\{\{(?:example:|case-gallery|manualhome|parameter-directory)',text(main_content)),
                rel+': unresolved handbook directive')
    require(used_downloads==set(downloads),'Example download manifest has missing or extra rows')
    all_html = set(identifiers)
    for lang in LANGUAGES:
        existing = {rel[len(lang)+1:] for rel in all_html if rel.startswith(lang+'/')}
        require(expected_pages<=existing,lang+': missing canonical pages '+str(sorted(expected_pages-existing)))
        require(set(output_rows[lang])==set(rd),lang+': reference-output topic inventory mismatch')
        search = read(site/lang/'search.json')
        indexed = {urlsplit(r['url']).path for r in search}
        require(indexed==expected_pages,lang+': wrong search coverage '+str(sorted(indexed ^ expected_pages)))
        function_hits = {r['title'] for r in search if r.get('kind')=='function'}
        require(exports<=function_hits,lang+': search omits exported functions '+str(sorted(exports-function_hits)))
        for row in search:
            remember(lang+'/search.json',row['url'])
        require({unquote(urlsplit(r['url']).fragment) for r in search if r.get('kind')=='parameter'} ==
                {'param-'+api+'-'+arg for api,arg in formals},lang+': parameter search results differ from public formals')
        counts['search_results_'+lang] = len(search)
    for page in expected_pages:
        require(body_anchors.get('en/'+page)==body_anchors.get('zh/'+page),page+': mismatched bilingual article anchors')
    pairs = {}
    for topic in rd:
        if not all((lang,topic) in references for lang in LANGUAGES):
            errors.append(topic+': missing native clipboard input')
            continue
        pairs[topic] = [references[(lang,topic)] for lang in LANGUAGES]
        en_sigs,zh_sigs = [signatures[(lang,topic)] for lang in LANGUAGES]
        require(len(en_sigs)==len(zh_sigs),topic+': different signature counts')
        for i,pair in enumerate(zip(en_sigs,zh_sigs)):
            pairs[topic+' signature '+str(i)] = pair
        require(figures[('en',topic)]==figures[('zh',topic)],topic+': bilingual example plots differ')
    scripts_for_comparison(pairs)
    version = pages['migration']['canonical_path']+'#api'
    for old in RETIRED:
        key = 'reference/'+old+'.html'
        require(key in migrations,key+': missing retired API version route')
        for prefix in ('','en/','zh/'):
            start = prefix+key
            end = (prefix or 'en/')+version
            for fragment in ('','#ref-examples'):
                cases.add((start+fragment,end))
            require(start in scripts,start+': missing version redirect')
            fallback = Document(site/start).root if (site/start).is_file() else ET.Element('body')
            targets = {destination(start,a.get('href','')) for a in fallback.iter('a')}
            langs = [prefix[:2]] if prefix else LANGUAGES
            for language in langs:
                require((language+'/'+pages['migration']['canonical_path'],'api') in targets,
                        start+': version page unavailable with JavaScript disabled ('+language+')')
    for lang in LANGUAGES:
        for old,mapping in migrations.items():
            source = lang+'/'+old
            if old.startswith('reference/') and Path(old).stem in RETIRED:
                continue
            for frag, dest in mapping.items():
                if not dest:
                    continue
                start = source+('#'+frag if frag else '')
                end = lang+'/'+dest
                if start==end:
                    continue
                remember(lang+'/index.html', dest)
                cases.add((start,end))
                require(source in scripts,source+': registered migration lacks a redirect')
    published = {r['path'] for r in read(ROOT/'tools/site-audit/baseline.json')['urls']}
    for name in published | set(read(WEB/'historical-pages.json')):
        require(name in all_html, name+': missing historical URL')
    for path,dest in routes.items():
        if path=='404.html':
            continue
        resolved = remember('index.html',dest)
        require(path in all_html and path in scripts,path+': missing legacy redirect')
        if resolved is None:
            continue
        name,fragment = resolved
        prefix,_,local = name.partition('/')
        if prefix in LANGUAGES and local in migrations:
            mapping = migrations[local]
            final = mapping.get(fragment, mapping.get('') if local in alias_pages else None)
            if final:
                dest = prefix+'/'+final
        cases.add((path,dest))
    for name in retired:
        require(name in scripts,name+': historical page must redirect')
        if name.startswith('reference/'):
            cases.add((name, 'en/'+pages['migration']['canonical_path']))
    for origin,url,name,anchor in links:
        dest = site/name
        require(dest.is_file(),origin+': missing local URL '+url)
        if anchor and name.endswith('.html'):
            require(anchor in identifiers.get(name,set()),origin+': missing fragment '+url)
    if not errors:
        run_redirects(scripts, sorted(cases))
    counts.update(handbook_pages=len(pages), native_topics=len(rd), handbook_inputs=len(records),
                  unique_examples=len(specs), handbook_plots=len(plots), local_urls=len(links), redirect_cases=len(cases),
                  native_figures_en=sum(len(v) for (lang,_),v in figures.items() if lang=='en'),
                  native_figures_zh=sum(len(v) for (lang,_),v in figures.items() if lang=='zh'))
    if args.report:
        args.report.write_text(json.dumps(dict(source_commit=manifest.get('source_commit'),counts=counts,errors=errors),ensure_ascii=False,indent=2)+'\n')
    if errors:
        raise SystemExit('\n'.join(dict.fromkeys(errors)))
    print('PASS:',json.dumps(counts,sort_keys=True))


if __name__=='__main__':
    main()
