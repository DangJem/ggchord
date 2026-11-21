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
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[2]
WEB = ROOT/'website'
LANGUAGES = ('en','zh')
VOID = set('area base br col embed hr img input link meta param source track wbr'.split())
BASE = 'https://local.invalid/ggchord/'
MARKER = re.compile(r'\{\{example:([a-z0-9-]+)\}\}')
PREVIEW = re.compile(r'\{\{preview:([a-z0-9-]+)\}\}')


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


def console_text(value):
    # Template rendering removes the padding printed after table columns.
    return '\n'.join(line.rstrip(' \t') for line in value.split('\n'))


def styled(elem, css, tag=None):
    return [n for n in elem.iter(tag) if css in n.get('class','').split()]


def check_article_directories(doc, main_content, rel, lang, local, navigation, pages, require):
    layouts = styled(doc, 'doc-layout')
    require(len(layouts)==1,rel+': expected a single article layout')
    if len(layouts)!=1:
        return
    layout = layouts[0]
    chapters = styled(layout, 'chapter-nav')
    outlines = styled(layout, 'page-nav')
    if local == 'index.html':
        require(not chapters and not outlines, rel+': homepage contains article directories')
        return
    expected_chapters = [group['id'] for group in navigation]
    require(len(chapters)==1, rel+': expected a single chapter directory')
    if len(chapters)!=1:
        return
    tree = chapters[0]
    groups = styled(tree, 'chapter-group')
    require([g.get('data-chapter') for g in groups]==expected_chapters, rel+': missing or reordered site chapters')
    current = [a for a in tree.iter('a') if a.get('aria-current')=='page']
    require(len(current)==1 and target(rel,current[0].get('href',''))==(rel,''),
            rel+': chapter directory must identify this page')
    group = next((g for g in groups if current and current[0] in g.iter('a')),None)
    expanded = [g for g in groups if any('open' in d.attrib for d in g.findall('./details'))]
    require(expanded==([group] if group is not None and group.find('./details') is not None else []),
            rel+': initial chapter expansion does not match the current page')
    expected_pages = {lang+'/'+page['canonical_path'] for page in pages.values()} | {rel}
    tree_links = [target(rel,a.get('href','')) for a in tree.iter('a')]
    require(set(tree_links)=={(page,'') for page in expected_pages} and len(tree_links)==len(expected_pages),
            rel+': chapter links differ from the handbook registry or repeat entries')

    excluded = {h for box in styled(main_content, 'lesson-extension', 'details') for h in box.iter()}
    levels = {'h2'} if local=='gallery.html' else {'h2','h3'}
    headings = [h for h in main_content.iter() if h.tag in levels and h.get('id') and h not in excluded]
    require(len(outlines)==int(bool(headings)),rel+': outline must be omitted when there are no article headings')
    if headings and outlines:
        links = list(outlines[0].iter('a'))
        section = subsection = 0
        numbers = []
        for heading in headings:
            if heading.tag=='h2' or not numbers:
                section += 1
                subsection = 0
            else:
                subsection += 1
            numbers.append(f'{section}.{subsection}' if subsection else str(section))
        prefixes = [a.find('./span[@class="rail-number"]') for a in links]
        require([a.get('href') for a in links]==['#'+h.get('id') for h in headings] and
                [text(prefix) for prefix in prefixes]==numbers and
                [text(a).removeprefix(n).strip() for a,n in zip(links,numbers)]==[text(h).strip() for h in headings],
                rel+': outline must follow article heading order and labels')
    require([node for node in layout if node in chapters+outlines+[main_content]]==chapters+outlines+[main_content],
            rel+': expected chapter directory, page outline and article in reading order')


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
    menu_sections = read(WEB/'manual-sections.json')
    nav_ui = read(WEB/'navigation-ui.json')
    groups = [g for g in navigation if g.get('navbar')]
    specs = read(WEB/'handbook-examples.json')
    home = read(WEB/'home.json')
    previews = {row['id']:row for row in read(WEB/'home-previews.json')}
    full_cases = {row['id']:row['page'] for row in read(WEB/'manual-figures.json')}
    example_cases = {p['canonical_path']:MARKER.findall((WEB/'en'/(p['slug']+'.md')).read_text())
                     for p in pages.values() if p['kind']=='case' and p['slug'].startswith('case-')}
    flex = read(WEB/'flexible-parameters.json')
    parameters = read(site/'parameter-directory.json')
    formals = {(r['api'],r['parameter']):r for r in parameters}
    manual = read(WEB/'reference-manual/structure.json')
    api = manual['topics']
    rd = {p.stem:re.findall(r'\\alias\{([^}]+)\}',p.read_text()) for p in (ROOT/'man').glob('*.Rd')}
    exports = set(re.findall(r'^export\(([^)]+)\)',(ROOT/'NAMESPACE').read_text(),re.M))
    records = read(site/'assets/examples.json')
    downloads = {(r['language'],r['page'],r['index']):r for r in records}
    rendered = read(site/'assets/handbook/rendered.json')
    output = read(site/'reference-output.json')
    output_rows = output['topics']
    migrations = read(WEB/'page-migrations.json')
    retired_api = {path for path in migrations if path.startswith('reference/')}
    routes = read(site/'legacy-routes.json')
    retired = set(read(site/'retired-routes.json'))
    identifiers, scripts, references, signatures = {}, {}, {}, {}
    links, cases, errors = set(), set(), []
    used_downloads = set()
    counts = Counter()
    figures = {}
    argument_terms = {}
    body_anchors = {}
    expected_pages = {p['canonical_path'] for p in pages.values()}
    canonical_topics = {'reference/'+slug+'.html' for slug in api}
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

    require(len(formals)==len(parameters) and all(r['api'] == r['topic'] and r['api'] in api for r in parameters),
            'Parameter directory names a duplicate or unknown exported argument')
    require(all((api,arg) in formals for api,args_here in flex.items() for arg in args_here),
            'Flexible argument missing from the runtime parameter directory')
    plots = {key for key,s in specs.items() if s['kind']=='plot'}
    png_dimensions = {}
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
            if len(header)==24:
                png_dimensions[name] = tuple(map(str,struct.unpack('>II', header[16:24])))
            dims = tuple(round(specs[name][d]*specs[name].get('dpi',144)) for d in ('width','height'))
            require(len(header)==24 and header[:8]==b'\x89PNG\r\n\x1a\n' and
                    all(abs(x-y)<=1 for x,y in zip(struct.unpack('>II', header[16:24]), dims)),
                    name+': PNG format/dimensions differ from plotting manifest')
    require(set(output_rows) == set(api), 'Reference output inventory differs from the published topics')
    require(output.get('status')=='passed', 'Native render reported failed example output')
    install_inputs = {name for name in manifest['source_sha256'] if name.startswith(('R/','data/','inst/','src/'))
                      or (name.startswith('man/') and name.endswith('.Rd'))}
    install_inputs.update(('DESCRIPTION','NAMESPACE','LICENSE','LICENSE.md'))
    with ZipFile(site/'assets/ggchord-source.zip') as archive:
        require(set(archive.namelist())=={'ggchord/'+name for name in install_inputs},
                'Matching package download omits inputs or includes extraneous files')
        for name in install_inputs:
            require(archive.read('ggchord/'+name)==(ROOT/name).read_bytes(),
                    'Download differs from the package used in the manual: '+name)
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
        require(not any(pre.find('.//a') is not None for pre in doc.iter('pre')), rel+': code contains a link')
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
        check_article_directories(doc,main_content,rel,lang,local,navigation,pages,require)
        body_anchors[rel] = main_ids
        require({'main','site-search','site-language'} <= identifiers[rel], rel+': missing accessibility, search or language entry')
        nav = styled(doc,'site-primary-nav')
        targets = [destination(rel,a.get('href','')) for a in styled(nav[0],'nav-link','a')] if nav else []
        require(targets==[(lang+'/'+n,'') for n in expected_nav],rel+': incorrect primary navigation')
        actual_menus = {n.get('id'):[destination(rel,a.get('href','')) for a in n.iter('a')]
                        for n in styled(doc,'site-nav-menu')}
        require(actual_menus=={key:[(lang+'/'+path,'') for path in values] for key,values in menus.items()},
                rel+': dropdown destinations differ from the navigation registry')
        for menu in styled(doc, 'site-nav-menu'):
            group_id = menu.get('id').removeprefix('site-menu-')
            sections = menu_sections.get(group_id, [])
            landing = next(g['landing'] for g in groups if g['id']==group_id)
            sections = [s for s in sections if any(key!=landing for key in s['pages'])]
            columns = menu.findall('./li')
            require(len(columns)==len(sections), rel+': incorrect menu group count for '+group_id)
            for column, section in zip(columns, sections):
                expected = [(lang+'/'+pages[k]['canonical_path'],'') for k in section['pages'] if k!=landing]
                require('menu-group' in column.get('class','').split() and
                        text(column.find('./span'))==section['title'][lang] and
                        [destination(rel,a.get('href','')) for a in column.findall('./ul/li/a')]==expected,
                        rel+': incorrect menu column '+section['title'][lang])
                for key, label in section.get('subheading_before',{}).items():
                    rows = [li for li in column.findall('./ul/li') if text(li.find('./span'))==label[lang]]
                    require(len(rows)==1 and rows[0].find('a') is not None and
                            destination(rel,rows[0].find('a').get('href',''))==(lang+'/'+pages[key]['canonical_path'],''),
                            rel+': missing subsection heading for '+key)
        utilities = styled(doc,'site-nav-tools')
        controls = utilities[0].findall('./li') if utilities else []
        require(len(controls)==4, rel+': incorrect toolbar control count')
        if len(controls)==4:
            search = controls[0].find('.//input[@id="site-search"][@type="search"][@aria-controls="site-results"]')
            github = controls[1].find('./a[@aria-label="GitHub"][@href]')
            appearance = controls[2].find('./button[@id="dropdown-lightswitch"][@aria-labelledby="appearance-label"]')
            language = controls[3].find('./details/summary')
            require(all(control is not None and control.get('tabindex') is None for control in
                        (search, github, appearance, language)), rel+': toolbar order or native focus is incorrect')
            require(search is not None and search.get('placeholder')==nav_ui['search_hint'][lang] and
                    text(controls[0].find('.//label[@for="site-search"]'))==nav_ui['search'][lang] and
                    controls[0].find('.//ul[@id="site-results"][@role="listbox"]') is not None,
                    rel+': inline search is incorrect')
            require(appearance is not None and text(appearance.find('./span[@id="appearance-label"]'))==nav_ui['appearance'][lang],
                    rel+': appearance label is incorrect')
            require(text(language).strip()=='🌐 '+nav_ui['language'][lang], rel+': visible language label is incorrect')
        other = 'zh' if lang=='en' else 'en'
        switch = doc.find('.//a[@id="site-language"]')
        require(switch is not None and destination(rel,switch.get('href',''))==(other+'/'+local,''),
                rel+': language switch does not lead to corresponding page')
        toggles = styled(doc,'site-nav-toggle','button')
        require({t.get('aria-controls') for t in toggles}==set(menus),rel+': missing menu control')
        for toggle in toggles:
            require(toggle.get('aria-controls') in identifiers[rel] and toggle.get('aria-expanded')=='false',
                    rel+': menu button has no matching controlled menu')
        if local == 'data.html':
            definitions = styled(main_content, 'field-definitions', 'div')
            require(len(definitions)==5 and all(
                'table-scroll' in block.get('class','').split() and
                len(block.findall('./table/thead/tr/th'))==3 for block in definitions),
                rel+': expected five independently sized field-definition tables')
            require(len(list(main_content.iter('table')))==len(definitions)+3 and
                    main_content.find('.//table/colgroup/col[@style]') is None,
                    rel+': data records or field descriptions still have inferred column widths')
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
            blocks = [b for b in styled(example,'code-disclosure','details') if b.get('data-code-kind')=='r']
            require(len(blocks)==1 and blocks[0].get('data-code-kind')=='r',rel+': expected one native R example')
            if not blocks:
                continue
            copied = clipboard(blocks[0])
            source = site/'assets/reference-examples'/(topic+'.R')
            require(source.is_file() and copied == source.read_text(), rel+': example differs from its download')
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
            row = output_rows.get(topic,{})
            console = styled(example,'data-output')
            require(imgs == ['../../assets/reference-examples/'+f for f in row['figures']] and
                    (bool(imgs) or bool(console)) and
                    [console_text(text(b.find('.//pre/code'))) for b in console] ==
                    ([console_text(row['console'])] if row['console'] else []),
                    rel+': example output differs from its R results')
            marked = set()
            sections_by_kind = {'sequence':'scalar','strand':'width','strand-offset':'strand',
                                'fixed-label':'scope','band':'scope'}
            label = '【支持灵活参数】' if lang=='zh' else '[Flexible parameter]'
            for term in main_content.iter('dt'):
                for name in styled(term,'argument-term'):
                    code = name.find('code')
                    arg = text(code)
                    argument_terms[(lang,topic,arg)] = {code.get('id')} if code is not None else set()
                    badges = styled(name,'flexible-parameter','a')
                    kind = flex.get(topic,{}).get(arg)
                    require(len(badges)==int(bool(kind)), f'{rel}: badge count for {arg}')
                    if badges and kind:
                        marked.add(arg)
                        badge = badges[0]
                        require(list(name)[-2:]==[code,badge] and text(badge)==label and
                                badge.get('href')=='../parameters.html#'+sections_by_kind[kind],
                                f'{rel}: wrong inline label or destination for {arg}')
            require(len(styled(main_content,'flexible-parameter'))==len(marked),rel+': stray flexible badge')
            require(marked==set(flex.get(topic,{})),f'{rel}: incorrect flexible parameter badges {marked}')
            if lang=='zh':
                require(bool(re.search('[\u3400-\u9fff]',text(main_content.find('.//h1')))) or
                        text(main_content.find('.//h1')).strip() == topic+'()', rel+': missing Han topic title')
            continue
        page = next((p for p in pages.values() if p['canonical_path']==local),None)
        if page is None:
            continue
        slug = page['slug']
        heading = ' '.join(text(main_content.find('.//h1')).split())
        expected_title = home['hero']['title'][lang] if slug=='index' else page['title'][lang]
        require(heading==expected_title,rel+': H1 differs from the registered title: '+heading)
        if slug=='reference':
            entries = {destination(rel,a.get('href','')) for group in styled(main_content,'reference-group')
                       for dl in group.iter('dl') for a in dl.iter('a')}
            require(entries=={(lang+'/'+key,'') for key in canonical_topics}, rel+': reference index mismatch')
            require('reference-filter' in identifiers[rel],rel+': no reference filter')
        if slug=='parameters':
            require(heading==('Flexible parameters' if lang=='en' else '灵活参数'),rel+': incorrect flexible-parameter guide heading')
        source = (WEB/lang/(slug+'.md')).read_text()
        expected = MARKER.findall(source)
        previews_here = PREVIEW.findall(source)
        visual_previews = styled(main_content,'handbook-preview','figure')
        preview_images = [destination(rel,img.get('src','')) for fig in visual_previews for img in fig.iter('img')]
        require(preview_images==[('assets/handbook/'+name+'.png','') for name in previews_here],
                rel+': introductory figures differ from their registered examples')
        for name,figure in zip(previews_here,visual_previews):
            spec = specs[name]
            caption = spec['caption'][lang]
            img = figure.find('.//img')
            link = figure.find('a')
            require(img is not None and img.get('alt')==caption and text(figure.find('figcaption')).strip()==caption
                    and link is not None and destination(rel,link.get('href',''))==
                    (lang+'/'+spec['preview_page']+'.html','example-'+name),
                    rel+': preview caption or code link differs for '+name)
        if slug=='installation':
            require(bool(previews_here) and not expected,rel+': the introduction should illustrate the steps without plotting code')
            require(not any(re.search(r'\b(ggchord|geom_\w+|coord_\w+)\s*\(',text(pre)) for pre in main_content.iter('pre')),
                    rel+': plotting code embedded in the introduction')
            require('install.packages("ggchord")' in text(main_content) and
                    any(destination(rel,a.get('href',''))==('assets/ggchord-source.zip','') for a in main_content.iter('a')),
                    rel+': missing CRAN or matching-source installation path')
        blocks = [n for n in main_content.iter() if n.get('data-example-kind')]
        dl_links = styled(main_content,'example-download','a')
        require(len(expected)==len(blocks)==len(dl_links),f'{rel}: marker/code/download counts {len(expected)}/{len(blocks)}/{len(dl_links)}')
        actual_figs = [destination(rel,img.get('src','')) for fig in styled(main_content,'handbook-output','figure') for img in fig.iter('img')]
        for fig in [*visual_previews, *styled(main_content,'handbook-output','figure')]:
            for img in fig.iter('img'):
                key = Path(img.get('src','')).stem
                require((img.get('width'),img.get('height'))==png_dimensions.get(key),rel+': HTML dimensions differ from PNG '+key)
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
            hero = home['hero']
            home_actions = [a for block in styled(main_content,'home-actions') for a in block.iter('a')]
            require([(a.get('href'),text(a)) for a in home_actions]==
                    [(a['url'],a['label'][lang]) for a in hero['actions']] and len(home_actions)==2,
                    rel+': homepage should have two introduction links')
            require([text(p).strip() for p in styled(main_content,'home-description','p')]==[hero['description'][lang]] and
                    len(styled(main_content,'home-logo','img'))==1,rel+': homepage introduction differs from its source')
            require({'choose','chord','circular','collinear','learn','quick-start','results','versions','modules'} <= main_ids,
                    rel+': missing a homepage bookmark')
            cards = styled(main_content,'home-module','section')
            modules = home['modules']
            require([card.get('data-module') for card in cards]==['multi-sequence','circular','dotplot','fusion'] and
                    [len(item['figures']) for item in modules]==[4,1,1,1],
                    rel+': expected four homepage groups and seven previews')
            for card,item in zip(cards,modules):
                pictures = styled(card,'home-preview','figure')
                keys = [Path(img.get('src','')).stem for preview in pictures for img in preview.iter('img')]
                require(keys==item['figures'],rel+': wrong flagship preview in '+item['slug'])
                anchor_id = 'home-'+item['slug']
                require(card.get('aria-labelledby')==anchor_id and
                        text(card.find(f'./h2[@id="{anchor_id}"]')).strip()==item['title'][lang],
                        rel+': missing or incorrect homepage group heading '+anchor_id)
                for key,preview in zip(keys,pictures):
                    case_links = styled(preview,'home-case','a')
                    links_to_figure = styled(preview,'figure-link','a')
                    require(len(case_links)==len(links_to_figure)==1,rel+': preview needs distinct case and full figure links')
                    if not (case_links and links_to_figure and key in previews):
                        continue
                    case = full_cases[previews[key]['full']]
                    ids = example_cases.get(case,[])
                    anchor = links_to_figure[0]
                    require(case_links[0].get('href')==anchor.get('data-code-url')==case and
                            text(case_links[0])==home['figure_labels']['case'][lang] and
                            text(anchor.find('span'))==home['figure_labels']['full'][lang] and len(ids)==1 and
                            destination(rel,anchor.get('href',''))==('assets/handbook/'+ids[0]+'.png',''),
                            rel+': homepage expanded image does not match case script for '+key)
                    title = home['figure_titles'][key][lang]
                    img = anchor.find('img')
                    require(img is not None and img.get('alt')==title and
                            text(preview.find('./figcaption/strong')).strip()==title,rel+': wrong homepage figure title '+key)
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
        require(not re.search(r'\{\{(?:example:|preview:|case-gallery|manualhome)',text(main_content)),
                rel+': unresolved handbook directive')
    require(used_downloads==set(downloads),'Example download manifest has missing or extra rows')
    all_html = set(identifiers)
    for lang in LANGUAGES:
        existing = {rel[len(lang)+1:] for rel in all_html if rel.startswith(lang+'/')}
        require(expected_pages<=existing,lang+': missing canonical pages '+str(sorted(expected_pages-existing)))
        for (function,arg),row in formals.items():
            require(row['anchor'] in argument_terms.get((lang,row['topic'],arg),set()),
                    f'{lang}: {function}({arg}) anchor does not belong to the documented argument')
        search = read(site/lang/'search.json')
        indexed = {urlsplit(r['url']).path for r in search}
        require(indexed==expected_pages,lang+': wrong search coverage '+str(sorted(indexed ^ expected_pages)))
        function_hits = [r['url'] for r in search if r.get('kind')=='function']
        require(sorted(function_hits)==sorted(canonical_topics),lang+': missing or duplicate function search results')
        parameter_hits = [r['url'] for r in search if r.get('kind')=='parameter']
        expected_hits = ['reference/'+r['topic']+'.html#'+r['anchor'] for r in parameters]
        require(sorted(parameter_hits)==sorted(expected_hits),lang+': parameter search destinations differ from the runtime formals')
        for row in search:
            remember(lang+'/search.json',row['url'])
        counts['search_results_'+lang] = len(search)
    for page in expected_pages:
        require(body_anchors.get('en/'+page)==body_anchors.get('zh/'+page),page+': mismatched bilingual article anchors')
    pairs = {}
    for topic in api:
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
    for key in retired_api:
        for prefix in ('','en/','zh/'):
            start = prefix+key
            for fragment in ('','#ref-examples'):
                cases.add((start+fragment,(prefix or 'en/') + migrations[key]['']))
            require(start in scripts,start+': missing documentation redirect')
            fallback = Document(site/start).root if (site/start).is_file() else ET.Element('body')
            targets = {destination(start,a.get('href','')) for a in fallback.iter('a')}
            for language in ([prefix[:2]] if prefix else LANGUAGES):
                require(target('index.html',language+'/'+migrations[key]['']) in targets,
                        start+': missing localized fallback for '+language)
    for lang in LANGUAGES:
        for old,mapping in migrations.items():
            source = lang+'/'+old
            if old in retired_api:
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
    removed_parameter_index = {'parameter-atlas.html', 'en/parameter-atlas.html', 'zh/parameter-atlas.html'}
    require(not removed_parameter_index & all_html, 'The removed parameter index is present in the build')
    for name in (published | set(read(WEB/'historical-pages.json'))) - removed_parameter_index:
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
    existing_reference_names = {alias for aliases in rd.values() for alias in aliases} | set(api) | {'index'}
    for name in read(WEB/'historical-pages.json'):
        if not name.startswith('reference/') or name in retired_api or Path(name).stem in existing_reference_names:
            continue
        fallback = Document(site/name).root if (site/name).is_file() else ET.Element('body')
        targets = {destination(name,a.get('href','')) for a in fallback.iter('a')}
        for lang in LANGUAGES:
            selected = [t for t in targets if t and t[0]==lang+'/'+pages['migration']['canonical_path']]
            require(len(selected)==1,name+': missing localized fallback to version selection for '+lang)
            if lang=='en' and selected:
                filename,fragment = selected[0]
                cases.add((name,filename+('#'+fragment if fragment else '')))
    for origin,url,name,anchor in links:
        dest = site/name
        require(dest.is_file(),origin+': missing local URL '+url)
        if anchor and name.endswith('.html'):
            require(anchor in identifiers.get(name,set()),origin+': missing fragment '+url)
    if not errors:
        run_redirects(scripts, sorted(cases))
    counts.update(handbook_pages=len(pages), api_topics=len(api), handbook_inputs=len(records),
                  unique_examples=len(specs), handbook_plots=len(plots), parameter_rows=len(formals),
                  ellipsis_rows=sum(arg=='...' for api,arg in formals), dotted_argument_rows=sum('.' in arg and arg!='...' for api,arg in formals), local_urls=len(links), redirect_cases=len(cases),
                  api_figures_en=sum(len(v) for (lang,_),v in figures.items() if lang=='en'),
                  api_figures_zh=sum(len(v) for (lang,_),v in figures.items() if lang=='zh'))
    if args.report:
        args.report.write_text(json.dumps(dict(source_commit=manifest.get('source_commit'),counts=counts,errors=errors),ensure_ascii=False,indent=2)+'\n')
    if errors:
        raise SystemExit('\n'.join(dict.fromkeys(errors)))
    print('PASS:',json.dumps(counts,sort_keys=True))


if __name__=='__main__':
    main()
