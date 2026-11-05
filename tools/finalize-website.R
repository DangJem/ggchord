# Attach paired-language navigation and local progressive enhancements to pkgdown pages.
args<-commandArgs(trailingOnly=TRUE)
publication<-identical(args,'publication')
version <- read.dcf('DESCRIPTION', fields = 'Version')[1, 1]
topic_files<-list.files('man','\\.Rd$',full.names=TRUE)
alias_topics<-unlist(lapply(topic_files,function(path) {
 rd<-tools::parse_Rd(path)
 aliases<-vapply(Filter(function(x) identical(attr(x,'Rd_tag'),'\\alias'),rd),
  function(x) paste(unlist(x),collapse=''),'')
 setNames(rep(sub('\\.Rd$','',basename(path)),length(aliases)),aliases)
}))
navigation <- jsonlite::fromJSON('website/navigation.json', simplifyVector=FALSE)
page_index <- jsonlite::fromJSON('website/pages.json', simplifyVector=FALSE)
section_numbers <- unlist(lapply(navigation, function(group)
  setNames(seq_along(group$pages), unlist(group$pages))), use.names = TRUE)
for(lang in c('en','zh')) {
 root<-file.path('.site',lang)
 for(path in list.files(root,'\\.html$',recursive=TRUE,full.names=TRUE)) {
  relative<-substring(path,nchar(root)+2L)
  deep<-length(strsplit(relative,'/',fixed=TRUE)[[1]])-1L
  base<-paste0(paste(rep('../',deep+1L),collapse=''),'assets/')
  partner<-paste0(paste(rep('../',deep+1L),collapse=''),if(lang=='en')'zh/' else 'en/',relative)
  doc<-xml2::read_html(path)
  xml2::xml_set_attr(xml2::xml_find_first(doc,'//html'),'lang',if(lang=='zh')'zh-Hans' else 'en')
  og<-xml2::xml_find_first(doc,'//meta[@property="og:image"]')
  if(!inherits(og,'xml_missing')) xml2::xml_set_attr(og,'content','https://dangjem.github.io/ggchord/assets/logo.png')
  # Native dropdown IDs derive from translated text; keep paired URLs/ARIA stable.
  for(menu in xml2::xml_find_all(doc,'//*[@data-bs-toggle="dropdown"]')) {
    title <- trimws(xml2::xml_text(menu))
    id <- if(title %in% c('绘图章节', 'Plot chapters')) 'site-visualizations'
      else if(title %in% c('专题', 'Topics')) 'site-topics' else NULL
    if(is.null(id)) next
    old_id <- xml2::xml_attr(menu,'id')
    xml2::xml_set_attr(menu,'id',id)
    labelled <- xml2::xml_find_all(doc,'//*[@aria-labelledby]')
    for(node in labelled) if(identical(xml2::xml_attr(node,'aria-labelledby'),old_id))
      xml2::xml_set_attr(node,'aria-labelledby',id)
  }
  main<-xml2::xml_find_first(doc,'//main')
  xml2::xml_set_attr(main,'tabindex','-1')
  slug <- sub('\\.html$','',basename(relative))
  xml2::xml_set_attr(xml2::xml_find_first(doc,'//body'),'data-page',if(startsWith(relative,'reference/')) 'reference' else slug)
  layout <- xml2::xml_parent(main)
  xml2::xml_set_attr(layout,'class',paste(xml2::xml_attr(layout,'class'),'doc-layout'))
  prefix <- paste(rep('../',deep),collapse='')
  current_page <- if(startsWith(relative, 'reference/')) 'reference' else slug
  page_link <- function(key, label = NULL) {
    page <- Filter(function(x) identical(x$slug, key), page_index)[[1]]
    target <- if(key == 'reference') 'reference/index.html' else paste0(key, '.html')
    if(is.null(label)) label <- paste0(section_numbers[[key]], '. ', page$title[[lang]])
    paste0('<a ', if(key == current_page) 'aria-current="page" ' else '',
      'href="', prefix, target, '">', htmltools::htmlEscape(label), '</a>')
  }
  chapter_links <- paste(vapply(navigation, function(group) {
    links <- paste(vapply(group$pages, page_link, ''), collapse = '')
    paste0('<details class="chapter-section"', if(current_page %in% group$pages || (current_page == 'index' && group$id == 'techniques')) ' open' else '',
      '><summary>', htmltools::htmlEscape(group$title[[lang]]), '</summary>', links, '</details>')
  }, ''), collapse = '')
  nav_fragment <- function(html) xml2::xml_find_first(xml2::read_html(paste0('<html><body>',html,'</body></html>')),'//body/*')
  xml2::xml_add_sibling(main,nav_fragment(paste0('<nav class="chapter-nav" aria-label="',
    if(lang=='zh')'章节与参考' else 'Chapters and reference','"><p class="navigation-title">',
    if(lang=='zh')'章节与参考' else 'Chapters and reference','</p>',chapter_links,'</nav>')),.where='before')
  xml2::xml_add_child(main,nav_fragment(paste0('<details class="chapter-drawer"><summary>',
    if(lang=='zh')'章节与参考' else 'Chapters and reference','</summary><nav>',chapter_links,'</nav></details>')),.where=0L)

  current_group <- Filter(function(group) current_page %in% group$pages, navigation)
  if(length(current_group) && current_page != 'index' && !startsWith(relative, 'reference/')) {
    group <- current_group[[1]]
    group_title <- htmltools::htmlEscape(group$title[[lang]])
    sequence <- group$id %in% c('chord', 'collinear', 'plasmid')
    stage <- group_title
    if(sequence) {
      step <- match(current_page, group$pages)
      total <- length(group$pages)
      stage <- paste(group_title, '·', if(lang == 'zh') sprintf('第 %d 节，共 %d 节', step, total)
        else sprintf('Section %d of %d', step, total))
      prev <- if(step > 1) page_link(group$pages[[step - 1]], if(lang == 'zh') '← 上一节' else '← Previous section')
      next_link <- if(step < total) page_link(group$pages[[step + 1]], if(lang == 'zh') '下一节 →' else 'Next section →')
      end <- page_link('about', if(lang == 'zh') '回到章节目录' else 'Back to chapter index')
      xml2::xml_add_child(main, nav_fragment(paste0('<nav class="lesson-pagination" aria-label="',
        if(lang == 'zh') '章节进度' else 'Chapter navigation', '">', prev,
        if(step == total) end else next_link, '</nav>')))
    }
    xml2::xml_add_sibling(xml2::xml_find_first(main, './/h1'), nav_fragment(
      paste0('<p class="lesson-stage">', stage, '</p>')), .where = 'before')
  }
  if (current_page == "index" || (length(current_group) && !startsWith(relative, "reference/"))) {
    chapter_number <- if (current_page == "index") integer() else section_numbers[[current_page]]
    add_number <- function(node, number) {
      xml2::xml_add_sibling(xml2::xml_contents(node)[[1]], nav_fragment(paste0('<span class="section-number">',
        number, ' </span>')), .where = 'before')
    }
    if (length(chapter_number)) add_number(xml2::xml_find_first(main, './/h1'), paste0(chapter_number, '.'))
    counts <- integer(3)
    for (heading in xml2::xml_find_all(main, './/*[self::h2 or self::h3 or self::h4][@id]')) {
      level <- as.integer(sub('h', '', xml2::xml_name(heading))) - 1L
      counts[level] <- counts[level] + 1L
      counts[seq_along(counts) > level] <- 0L
      add_number(heading, paste(c(chapter_number, counts[seq_len(level)]), collapse = '.'))
    }
  }
  if(lang=='zh') {
   # Finalize native interface prose independently of runner gettext locales.
   ui <- c('Skip to contents'='跳到正文','Site navigation'='网站导航',
     'Toggle navigation'='展开导航','Light'='浅色','Dark'='深色','Auto'='自动',
     'On this page'='本页目录','Usage'='用法','Arguments'='参数',
     'Value'='返回值','Details'='详细说明','Examples'='示例','See also'='另见',
     'Format'='格式','Source'='来源','References'='参考资料','Author'='作者',
     'Note'='说明','Notes'='说明','Package index'='函数文档','包索引'='函数文档',
     'Authors and Citation'='作者与引用','Citation'='引用','Contents'='目录')
   for(node in xml2::xml_find_all(doc,'//text()[not(ancestor::pre or ancestor::script or ancestor::style)]')) {
    value <- trimws(xml2::xml_text(node))
    if(value %in% names(ui)) xml2::xml_set_text(node,ui[[value]])
    else if(grepl('^Developed by ',value)) xml2::xml_set_text(node,sub('Developed by ','开发者：',value,fixed=TRUE))
    else if(grepl('^Site built with',value)) xml2::xml_set_text(node,sub('Site built with','网站构建工具：',value,fixed=TRUE))
   }
   for(node in xml2::xml_find_all(doc,'//*[@aria-label or @title]')) for(attr in c('aria-label','title')) {
    value <- xml2::xml_attr(node,attr)
    if(!is.na(value) && value %in% names(ui)) xml2::xml_set_attr(node,attr,ui[[value]])
   }
   translations<-c('值'='返回值','范例'='示例','该页面'='本页目录','细节'='详细说明','相关话题'='另见')
   for(h in xml2::xml_find_all(doc,'//h2')) {
    children<-xml2::xml_contents(h)
    for(t in children) if(xml2::xml_type(t)=='text') {
     text<-trimws(xml2::xml_text(t))
     if(text %in% names(translations)) xml2::xml_set_text(t,translations[[text]])
    }
   }
  }
  # Native reference TOCs are populated later by JavaScript. Build the compact
  # outline from actual headings, rather than cloning their initially empty nav.
  outline_headings <- xml2::xml_find_all(main,
    './/*[self::h2 or self::h3][@id][not(ancestor::details[@class="lesson-extension"])]')
  if(length(outline_headings)) {
    outline_links <- paste(vapply(outline_headings,function(h) paste0('<li><a href="#',
      htmltools::htmlEscape(xml2::xml_attr(h,'id'),attribute=TRUE),'">',
      htmltools::htmlEscape(xml2::xml_text(h)),'</a></li>'),''),collapse='')
    toc_fragment <- nav_fragment(paste0('<details class="page-outline"><summary>',
      if(lang=='zh')'本页目录' else 'On this page','</summary><ul>',outline_links,'</ul></details>'))
    xml2::xml_add_sibling(xml2::xml_find_first(main,'.//h1'),toc_fragment,.where='after')
  }

  native_sidebar <- xml2::xml_find_all(layout,'./aside')
  xml2::xml_remove(native_sidebar)
  if(length(outline_headings)) {
    page_nav <- nav_fragment(paste0('<aside class="page-nav" aria-label="',
      if(lang=='zh')'本页目录' else 'On this page','"><p class="navigation-title">',
      if(lang=='zh')'本页目录' else 'On this page','</p><nav><ul>',outline_links,'</ul></nav></aside>'))
    xml2::xml_add_sibling(main,page_nav,.where='before')
  }

  for(logo in xml2::xml_find_all(main,'.//img[@src]')) {
   if(grepl('(^|/)logo\\.(svg|png)$',xml2::xml_attr(logo,'src'))) {
    xml2::xml_set_attr(logo,'src',paste0(base,'logo.webp'))
    xml2::xml_set_attr(logo,'alt','ggchord')
   }
  }
  head<-xml2::xml_find_first(doc,'//head')
  fragment<-function(text) xml2::xml_find_first(xml2::read_html(paste0('<html><body>',text,'</body></html>')),'//body/*')
  xml2::xml_add_child(head,fragment(paste0('<link rel="alternate" hreflang="',if(lang=='en')'zh-Hans' else 'en','" href="',partner,'">')))
  xml2::xml_add_child(head,fragment(paste0('<link rel="icon" href="',base,'favicon.png">')))
  # Add style and script separately so neither can be lost by first-node extraction.
  xml2::xml_add_child(head,fragment(paste0('<link rel="stylesheet" href="',base,'site.css?v=',unname(tools::md5sum('website/assets/site.css')),'">')))
  xml2::xml_add_child(head,fragment(paste0('<script src="',base,'site.js?v=',unname(tools::md5sum('website/assets/site.js')),'" defer></script>')))
  toggle<-xml2::xml_find_first(doc,'//button[contains(@class,"navbar-toggler")]')
  language<-fragment(paste0('<details class="language-control"><summary><span aria-hidden="true">🌐</span> ',
   if(lang=='zh')'现代汉语' else 'English',' <span aria-hidden="true">▾</span></summary>',
   '<nav aria-label="',if(lang=='zh')'选择语言' else 'Choose language','">',
   '<a ',if(lang=='en')'aria-current="true" ' else '', 'href="',if(lang=='en') '' else partner,
   '" lang="en" data-language="en">English</a>',
   '<a ',if(lang=='zh')'aria-current="true" ' else '', 'href="',if(lang=='zh') '' else partner,
   '" lang="zh-Hans" data-language="zh" ',if(lang=='en')'id="site-language"' else '',
   '>现代汉语</a></nav></details>'))
  # Keep one stable identifier for the other-language link on every page.
  if(lang=='zh') xml2::xml_set_attr(xml2::xml_find_first(language,'.//a[@lang="en"]'),'id','site-language')
  xml2::xml_add_sibling(toggle,language,.where='before')
  target<-xml2::xml_find_first(doc,'(//div[@id="navbar"]/ul[contains(@class,"navbar-nav")])[last()]')
  label<-if(lang=='zh')'搜索文档' else 'Search documentation'
  search<-fragment(paste0('<li class="nav-item local-search"><label class="visually-hidden" for="site-search">',label,
   '</label><input id="site-search" type="search" placeholder="',label,'" autocomplete="off" class="form-control">',
   '<div id="site-results" class="search-results" aria-live="polite"></div></li>'))
  xml2::xml_add_child(target,search)
  body<-xml2::xml_find_first(doc,'//body')
  xml2::xml_set_attr(body,'data-lang',lang)
  xml2::xml_set_attr(body,'data-search',paste0(paste(rep('../',deep),collapse=''),'search.json'))
  notice <- paste(version, if(lang == 'zh') '文档' else 'documentation')
  if(!publication) notice<-paste(notice,if(lang=='zh')'· 本地预览，尚未部署' else '· local preview, not deployed')
  xml2::xml_add_child(main,fragment(paste0('<div class="candidate-note">',notice,'</div>')),.where=0L)
  if(startsWith(relative,'reference/')) {
   info<-if(lang=='zh') '本页完整说明已译为现代汉语；函数、参数、代码和原始数据标识保留原样。' else 'Complete canonical API text; task articles explain the workflows.'
   link<-paste0(paste(rep('../',deep),collapse=''),'about.html')
   xml2::xml_add_child(main,fragment(paste0('<p class="reference-language-note">',info,' <a href="',link,'">',if(lang=='zh')'使用教程' else 'Tutorials','</a></p>')),.where=1L)
  }
  # Native locale messages add a translated S3 preamble to usage code. Keep
  # this generated comment canonical too, while translating surrounding prose.
  if(lang=='zh') for(comment in xml2::xml_find_all(main,'.//pre//span[@class="co"]')) {
   value <- xml2::xml_text(comment)
   if(grepl('^# .+ 类的 S3 方法$',value)) {
    class_name <- sub(' 类的 S3 方法$','',sub('^# ','',value))
    xml2::xml_set_text(comment,paste0("# S3 method for class '",class_name,"'"))
   }
  }
  for(pre in xml2::xml_find_all(main,'.//pre')) {
   parent<-xml2::xml_parent(pre)
   if(xml2::xml_name(parent)=='details')next
   details<-fragment(paste0('<details class="code-disclosure" open><summary>',if(lang=='zh')'代码' else 'Code','</summary>',as.character(pre),'</details>'))
   xml2::xml_add_sibling(pre,details,.where='before')
   xml2::xml_remove(pre)
  }
  # pkgdown rewrites Markdown links to HTML; these are downloadable source files.
  for(a in xml2::xml_find_all(main,'.//a[@href]')) {
   href<-xml2::xml_attr(a,'href')
   if(grepl('^https://github.com/DangJem/ggchord/blob/HEAD/R/',href))
    xml2::xml_set_attr(a,'href',paste0(base,'source/',basename(href)))
   if(grepl('/assets/(NEWS|ASSET-SOURCES|MIGRATION|MIGRATION-Hans)\\.html$',href))
    xml2::xml_set_attr(a,'href',sub('\\.html$','.md',href))
   href <- sub('^https://dangjem.github.io/ggchord/reference/', 'https://rdrr.io/pkg/ggchord/man/', href)
   if(grepl('^https?://rdrr\\.io/pkg/ggchord/man/[^#]+\\.html',href)) {
    name<-utils::URLdecode(sub('\\.html.*$','',sub('^.*/man/','',href)))
    topic<-if(name %in% names(alias_topics)) alias_topics[[name]] else
     if(name %in% sub('\\.Rd$','',basename(topic_files))) name else NA_character_
    if(is.na(topic)) stop('Unresolved current-package autolink: ',href)
    fragment<-if(grepl('#',href,fixed=TRUE)) sub('^[^#]+','',href) else ''
    xml2::xml_set_attr(a,'href',paste0(paste(rep('../',deep),collapse=''),'reference/',topic,'.html',fragment))
   }
  }

  # Translated native headings otherwise receive punctuation-only IDs. Match
  # paired heading positions to the finalized English page's stable URLs.
  if(lang=='zh') {
   english<-xml2::read_html(file.path('.site/en',relative))
   headings<-xml2::xml_find_all(main,'.//*[self::h2 or self::h3 or self::h4][@id]')
   counterparts<-xml2::xml_find_all(english,'//main//*[self::h2 or self::h3 or self::h4][@id]')
   if(length(headings)!=length(counterparts)) stop('Unmatched paired headings: ',relative)
   for(i in seq_along(headings)) {
    old<-xml2::xml_attr(headings[[i]],'id');fresh<-xml2::xml_attr(counterparts[[i]],'id')
    if(old==fresh)next
    links<-xml2::xml_find_all(doc,'//a[@href]')
    links<-links[xml2::xml_attr(links,'href')==paste0('#',old)]
    xml2::xml_set_attr(links,'href',paste0('#',fresh))
    xml2::xml_set_attr(headings[[i]],'id',fresh)
   }
  }
  # Native argument anchors can collide for names differing only by digits (e1/e2).
  # Preserve the first URL and make each following element's own anchor unique.
  nodes<-xml2::xml_find_all(doc,'//*[@id]')
  occupied<-xml2::xml_attr(nodes,'id');seen<-character()
  for(node in nodes) {
   id<-xml2::xml_attr(node,'id')
   if(id %in% seen) {
    i<-2L;fresh<-paste0(id,'-',i)
    while(fresh %in% occupied) {i<-i+1L;fresh<-paste0(id,'-',i)}
    xml2::xml_set_attr(node,'id',fresh)
    children<-xml2::xml_find_all(node,'.//a[@href]')
    children<-children[xml2::xml_attr(children,'href')==paste0('#',id)]
    xml2::xml_set_attr(children,'href',paste0('#',fresh))
    occupied<-c(occupied,fresh)
   }
   seen<-c(seen,id)
  }
  # API code is not executed, so the native search configuration is replaced by
  # a deterministic local full-text index with Han substring support.
  xml2::write_html(doc,path)
 }
}
cat('Paired navigation and local search attached to all documentation pages.\n')
