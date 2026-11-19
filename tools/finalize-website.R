# Attach paired-language navigation and local progressive enhancements to pkgdown pages.
args<-commandArgs(trailingOnly=TRUE)
publication<-identical(args,'publication')
`%||%` <- function(x,y) if(is.null(x) || length(x)==0L || is.na(x)) y else x
version <- read.dcf('DESCRIPTION', fields = 'Version')[1, 1]
topic_files<-list.files('man','\\.Rd$',full.names=TRUE)
alias_topics<-unlist(lapply(topic_files,function(path) {
 rd<-tools::parse_Rd(path)
 aliases<-vapply(Filter(function(x) identical(attr(x,'Rd_tag'),'\\alias'),rd),
  function(x) paste(unlist(x),collapse=''),'')
 setNames(rep(sub('\\.Rd$','',basename(path)),length(aliases)),aliases)
}))
source('tools/website/reference.R')
navigation <- jsonlite::fromJSON('website/navigation.json', simplifyVector=FALSE)
page_index <- jsonlite::fromJSON('website/pages.json', simplifyVector=FALSE)
page_migrations <- jsonlite::fromJSON('website/page-migrations.json', simplifyVector=FALSE)
page_record <- function(key) {
  found <- Filter(function(x) identical(x$slug,key),page_index)
  if(length(found)) found[[1]] else NULL
}
manual_sections <- jsonlite::fromJSON('website/manual-sections.json', simplifyVector=FALSE)
navigation <- setNames(navigation, vapply(navigation, function(group) group$id, ''))
for(lang in c('en','zh')) {
 root<-file.path('.site',lang)
 for(path in list.files(root,'\\.html$',recursive=TRUE,full.names=TRUE)) {
  relative<-substring(path,nchar(root)+2L)
  deep<-length(strsplit(relative,'/',fixed=TRUE)[[1]])-1L
  base<-paste0(paste(rep('../',deep+1L),collapse=''),'assets/')
  partner<-paste0(paste(rep('../',deep+1L),collapse=''),if(lang=='en')'zh/' else 'en/',relative)
  doc<-xml2::read_html(path)
  if(relative=='404.html') {
    locale_url <- paste0('https://dangjem.github.io/ggchord/',lang,'/')
    for(attribute in c('href','src')) {
      nodes <- xml2::xml_find_all(doc,paste0('//*[@',attribute,']'))
      links <- xml2::xml_attr(nodes,attribute)
      local <- startsWith(links,locale_url)
      xml2::xml_set_attr(nodes[local],attribute,substring(links[local],nchar(locale_url)+1L))
    }
  }
  xml2::xml_set_attr(xml2::xml_find_first(doc,'//html'),'lang',if(lang=='zh')'zh-Hans' else 'en')
  og<-xml2::xml_find_first(doc,'//meta[@property="og:image"]')
  if(!inherits(og,'xml_missing')) xml2::xml_set_attr(og,'content','https://dangjem.github.io/ggchord/assets/logo.png')
  for(navbar in xml2::xml_find_all(doc,'//nav[contains(@class,"navbar-expand-lg")]'))
    xml2::xml_set_attr(navbar,'class',sub('navbar-expand-lg','navbar-expand-xl',xml2::xml_attr(navbar,'class'),fixed=TRUE))
  main<-xml2::xml_find_first(doc,'//main')
  xml2::xml_set_attr(main,'tabindex','-1')
  slug <- sub('\\.html$','',basename(relative))
  if(startsWith(relative, 'reference/')) add_reference_navigation(doc, main, slug, lang)
  if(slug=='LICENSE') {
    label <- if(lang=='zh') 'MIT 许可' else 'MIT License'
    xml2::xml_set_text(xml2::xml_find_first(main,'.//h1'),label)
    xml2::xml_set_text(xml2::xml_find_first(doc,'//title'),paste(label,'• ggchord'))
  }
  if(slug=='404' && lang=='zh') for(node in xml2::xml_find_all(main,'.//p'))
    xml2::xml_set_text(node,'未找到页面。请从首页、绘图指南或搜索找回所需内容。')
  xml2::xml_set_attr(xml2::xml_find_first(doc,'//body'),'data-page',if(identical(relative,'reference/index.html')) 'reference' else slug)
  if(slug=='index' && length(xml2::xml_find_all(main,'.//h1'))>1L) xml2::xml_remove(xml2::xml_find_first(main,'./h1'))
  layout <- xml2::xml_parent(main)
  xml2::xml_set_attr(layout,'class',paste(xml2::xml_attr(layout,'class'),'doc-layout'))
  prefix <- paste(rep('../',deep),collapse='')
  current_page <- if(startsWith(relative, 'reference/')) 'reference' else slug
  record <- page_record(current_page)
  is_known_page <- !is.null(record)
  group <- if(is_known_page) navigation[[record$module]] else navigation$information
  module <- group$id
  primary <- if(is.null(group$parent)) group else navigation[[group$parent]]
  xml2::xml_set_attr(xml2::xml_find_first(doc,'//body'),'data-module',module)
  nav_fragment <- function(html) xml2::xml_find_first(xml2::read_html(paste0('<html><body>',html,'</body></html>')),'//body/*')
  escape <- htmltools::htmlEscape
  page_link <- function(key, label = NULL, class = NULL) {
    page <- page_record(key)
    if(is.null(label)) label <- page$title[[lang]]
    paste0('<a ',
      if(!is.null(class)) paste0('class="',class,'" '),
      if(identical(relative,page$canonical_path)) 'aria-current="page" ',
      'href="',prefix,escape(page$canonical_path,attribute=TRUE),'">',escape(label),'</a>')
  }

  # Keep each destination directly clickable and use a separate menu button.
  primary_groups <- Filter(function(x) !is.null(x$navbar), navigation)
  nav_items <- c(paste0('<li class="site-nav-item site-nav-home">',
    page_link('index',if(lang=='zh') '首页' else 'Home','nav-link'),'</li>'),
    vapply(primary_groups,function(g) {
      descendants <- Filter(function(x) identical(x$parent,g$id),navigation)
      keys <- unique(c(unlist(g$pages),unlist(lapply(descendants,function(x) c(x$landing,unlist(x$pages))))))
      keys <- Filter(function(key) !identical(key,g$landing) && !is.null(page_record(key)),keys)
      link <- page_link(g$landing,g$title[[lang]],'nav-link')
      if(identical(primary$id,g$id) && current_page!='index' && !identical(relative,page_record(g$landing)$canonical_path))
        link <- sub('<a ', '<a aria-current="location" ',link,fixed=TRUE)
      menu <- if(length(keys)) paste0('<button type="button" class="site-nav-toggle" aria-expanded="false" aria-controls="site-menu-',
        escape(g$id,attribute=TRUE),'" aria-label="',escape(paste(if(lang=='zh') '展开' else 'Open',g$title[[lang]]),attribute=TRUE),'">▾</button>',
        '<ul class="site-nav-menu" id="site-menu-',escape(g$id,attribute=TRUE),'">',
        paste(vapply(keys,function(key) paste0('<li>',page_link(key),'</li>'),''),collapse=''),'</ul>') else ''
      paste0('<li class="site-nav-item site-nav-',escape(g$navbar,attribute=TRUE),'">',link,menu,'</li>')
    },''))
  primary_list <- xml2::xml_find_first(doc,'//div[@id="navbar"]/ul[contains(@class,"navbar-nav")]')
  xml2::xml_set_attr(primary_list,'class','navbar-nav site-primary-nav')
  xml2::xml_remove(xml2::xml_children(primary_list))
  for(item in nav_items) xml2::xml_add_child(primary_list,nav_fragment(item))

  sections <- manual_sections[[module]]
  if(length(sections)) {
    topic_label <- if(lang=='zh') '相关主题' else 'Related topics'
    links <- vapply(sections,function(section) {
      items <- vapply(unlist(section$pages),function(key) paste0('<li>',page_link(key),'</li>'),'')
      paste0('<div class="manual-section"><span class="manual-section-title">',
        escape(section$title[[lang]]),'</span><ul class="manual-section-pages">',
        paste(items,collapse=''),'</ul></div>')
    },'')
    topics <- paste0(page_link(group$landing,group$title[[lang]],'module-home'),
      paste(links,collapse=''),
      if(!is.null(group$parent)) page_link(primary$landing,primary$title[[lang]],'topic-overview'))
    xml2::xml_add_sibling(main,nav_fragment(paste0('<nav class="chapter-nav" aria-label="',
      topic_label,'"><p class="navigation-title">',topic_label,'</p>',topics,'</nav>')),.where='before')
    xml2::xml_add_child(main,nav_fragment(paste0('<details class="chapter-drawer"><summary>',
      topic_label,'</summary><nav aria-label="',topic_label,'">',topics,'</nav></details>')),.where=0L)
  }

  if(current_page!='index') {
    heading <- xml2::xml_find_first(main,'.//h1')
    ancestors <- unique(c(primary$landing,group$landing))
    ancestors <- Filter(function(key) !identical(page_record(key)$canonical_path,relative), ancestors)
    crumbs <- c(page_link('index',if(lang=='zh') '首页' else 'Home'),
      vapply(ancestors,page_link,''),
      paste0('<span aria-current="page">',escape(xml2::xml_text(heading)),'</span>'))
    crumb <- paste0('<nav class="breadcrumb" aria-label="',if(lang=='zh') '位置导航' else 'Breadcrumb',
      '">',paste(crumbs,collapse='<span aria-hidden="true">/</span>'),'</nav>')
    xml2::xml_add_sibling(heading,nav_fragment(crumb),.where='before')
  }
  footer <- xml2::xml_find_first(doc,'//footer')
  if(!inherits(footer,'xml_missing')) {
    project <- page_link('about',class='project-link')
    xml2::xml_add_child(footer,nav_fragment(paste0('<p>',project,'</p>')))
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
  outline_headings <- xml2::xml_find_all(main, if(slug=='gallery')
    './/h2[@id][not(ancestor::details[@class="lesson-extension"])]' else
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
  brand <- xml2::xml_find_first(doc,'//a[contains(@class,"navbar-brand")]')
  xml2::xml_add_child(brand,nav_fragment(paste0('<img src="',base,'logo.webp" width="36" height="36" alt="" aria-hidden="true">')),.where=0L)
  head<-xml2::xml_find_first(doc,'//head')
  fragment<-function(text) xml2::xml_find_first(xml2::read_html(paste0('<html><body>',text,'</body></html>')),'//body/*')
  xml2::xml_add_child(head,fragment(paste0('<link rel="alternate" hreflang="',if(lang=='en')'zh-Hans' else 'en','" href="',partner,'">')))
  xml2::xml_add_child(head,fragment(paste0('<link rel="icon" href="',base,'favicon.png">')))
  # Add style and script separately so neither can be lost by first-node extraction.
  xml2::xml_add_child(head,fragment(paste0('<link rel="stylesheet" href="',base,'site.css?v=',unname(tools::md5sum('website/assets/site.css')),'">')))
  xml2::xml_add_child(head,fragment(paste0('<script src="',base,'site.js?v=',unname(tools::md5sum('website/assets/site.js')),'" defer></script>')))
  no_script <- xml2::xml_add_child(head,'noscript')
  xml2::xml_add_child(no_script,'style',paste0(
    '.navbar .navbar-toggler,.navbar .local-search,.navbar #dropdown-lightswitch{display:none}',
    '.site-nav-item:hover .site-nav-menu,.site-nav-item:focus-within .site-nav-menu{display:block}',
    '@media(max-width:1199.98px){',
    'body{padding-top:0}.navbar.fixed-top{position:static}',
    '.navbar .navbar-collapse.collapse{display:flex}',
    '.navbar .navbar-collapse>.navbar-nav:first-child{flex-direction:row;flex-wrap:wrap;gap:0 .35rem}',
    '}'))
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
  label<-if(lang=='zh')'搜索说明' else 'Search documentation'
  search<-fragment(paste0('<li class="nav-item local-search"><button type="button" class="search-trigger" id="site-search-trigger" aria-haspopup="dialog"><span aria-hidden="true">⌕</span><span class="search-trigger-label">',
   if(lang=='zh')'搜索文档' else 'Search docs','</span><kbd aria-hidden="true">⌘ K</kbd></button></li>'))
  xml2::xml_add_child(target,search)
  body<-xml2::xml_find_first(doc,'//body')
  xml2::xml_set_attr(body,'data-kind',if(startsWith(relative,'reference/')) 'api' else if(!is.null(record)) record$kind else 'utility')
  metadata <- jsonlite::fromJSON('.site/build-metadata.json')
  xml2::xml_remove(xml2::xml_find_all(doc,'//nav//small[contains(@class,"nav-text")]'))
  xml2::xml_add_child(body,fragment(paste0('<dialog class="search-dialog" id="site-search-dialog" aria-labelledby="site-search-title"><div class="search-heading"><h2 id="site-search-title">',label,
   '</h2><button type="button" class="search-close" aria-label="',if(lang=='zh')'关闭搜索' else 'Close search','">×</button></div>',
   '<label class="search-field" for="site-search"><span class="visually-hidden">',label,'</span><input id="site-search" type="search" role="combobox" aria-autocomplete="list" aria-controls="site-results" aria-expanded="false" autocomplete="off" placeholder="',
   if(lang=='zh')'例如：覆盖度、geom_feature、图例' else 'Examples: coverage, geom_feature, legends','"></label><div class="search-filters" aria-label="',
   if(lang=='zh')'筛选结果类型' else 'Filter result types','"></div><p class="search-count" id="search-status" role="status"></p><ul id="site-results" class="search-results" role="listbox" aria-label="',
   if(lang=='zh')'搜索结果' else 'Search results','"></ul><div class="search-help">',
   if(lang=='zh')'↑ ↓ 选择 · Enter 打开 · Esc 关闭；也可按任务查找。' else '↑ ↓ select · Enter open · Esc close; or browse by task.',
   ' <a href="',prefix,escape(page_record(navigation$plots$landing)$canonical_path,attribute=TRUE),'">',if(lang=='zh')'绘图指南' else 'Plotting guide','</a></div></dialog>')))
  xml2::xml_set_attr(body,'data-lang',lang)
  xml2::xml_set_attr(body,'data-search',paste0(paste(rep('../',deep),collapse=''),'search.json'))
  # Native locale messages add a translated S3 preamble to usage code. Keep
  # this generated comment canonical too, while translating surrounding prose.
  if(lang=='zh') for(comment in xml2::xml_find_all(main,'.//pre//span[@class="co"]')) {
   value <- xml2::xml_text(comment)
   if(grepl('^# .+ 类的 S3 方法$',value)) {
    class_name <- sub(' 类的 S3 方法$','',sub('^# ','',value))
    xml2::xml_set_text(comment,paste0("# S3 method for class '",class_name,"'"))
   }
  }
  examples <- xml2::xml_find_all(main, './/*[@data-example-kind]')
  for (i in seq_along(examples)) {
    name <- sprintf('%s-%d-%s.R', current_page, i, lang)
    if (!file.exists(file.path('.site/assets/examples', name))) stop('Missing example script: ',name)
    xml2::xml_add_sibling(examples[[i]], fragment(paste0('<a class="example-download" download href="',
      base,'examples/', name,'">',if(lang=='zh')'下载这段 R 代码' else 'Download this R example','</a>')), .where='after')
  }
  for(pre in xml2::xml_find_all(main,'.//pre')) {
   parent<-xml2::xml_parent(pre)
   if(xml2::xml_name(parent)=='details')next
   is_r <- grepl('sourceCode|r-output|usage',paste(xml2::xml_attr(pre,'class'),xml2::xml_attr(parent,'class')))
   usage <- length(xml2::xml_find_all(pre, "ancestor::div[h2[@id='ref-usage']]")) > 0L
   code_kind <- if (usage) 'signature' else if (is_r) 'r' else 'output'
   label <- switch(code_kind,
     signature = if(lang=='zh') '调用格式' else 'Usage syntax',
     r = if(lang=='zh') 'R · 示例代码' else 'R · Example code',
     output = if(lang=='zh') '输出 / 文本' else 'Output / text')
   details<-fragment(paste0('<details class="code-disclosure" data-code-kind="',code_kind,'" open><summary>',label,'</summary>',as.character(pre),'</details>'))
   xml2::xml_add_sibling(pre,details,.where='before')
   xml2::xml_remove(pre)
  }
  for(table in xml2::xml_find_all(main,'.//table')) {
    container <- xml2::xml_parent(table)
    if(grepl('data-output',xml2::xml_attr(container,'class') %||% '')) next
    if(length(xml2::xml_find_all(table,'(.//tr)[1]/*')) >= 4L)
      xml2::xml_set_attr(table,'class',trimws(paste(xml2::xml_attr(table,'class') %||% '', 'wide-table')))
    wrapper <- fragment(paste0('<div class="table-scroll" role="region" tabindex="0" aria-label="',if(lang=='zh')'可横向滚动的表格' else 'Horizontally scrollable table','"></div>'))
    xml2::xml_add_sibling(table,wrapper,.where='before')
    xml2::xml_add_child(xml2::xml_find_first(table,'preceding-sibling::*[1]'),table)
    xml2::xml_remove(table)
  }
  # Numeric alignment belongs to columns, including the header. Mixed
  # default-value columns remain left aligned.
  for(table in xml2::xml_find_all(main,'.//table')) {
    rows <- xml2::xml_find_all(table,'.//tbody/tr')
    if(!length(rows)) next
    n <- max(vapply(rows,function(row) length(xml2::xml_find_all(row,'./th|./td')),integer(1)))
    for(i in seq_len(n)) {
      cells <- xml2::xml_find_all(table,paste0('.//tbody/tr/*[self::td or self::th][',i,']'))
      values <- trimws(xml2::xml_text(cells)); values <- values[nzchar(values)]
      if(!length(values) || !all(grepl('^[+-]?([0-9]+([.,][0-9]+)?|\\.[0-9]+)([eE][+-]?[0-9]+)?%?$',values))) next
      headers <- xml2::xml_find_all(table,paste0('.//thead/tr/*[self::th or self::td][',i,']'))
      for(cell in c(headers,cells)) xml2::xml_set_attr(cell,'class',trimws(paste(xml2::xml_attr(cell,'class') %||% '', 'numeric')))
    }
  }
  for(output in xml2::xml_find_all(main,'.//*[contains(@class,"data-output")]')) {
    xml2::xml_set_attr(output,'tabindex','0')
    xml2::xml_set_attr(output,'role','region')
    xml2::xml_set_attr(output,'aria-label',if(lang=='zh')'R 数据输出，可横向滚动' else 'R data output, horizontally scrollable')
  }
  if(current_page=='index') xml2::xml_set_attr(xml2::xml_find_all(main,'.//*[contains(@class,"home-results")]//img'),
    'sizes','(max-width: 899px) calc(100vw - 48px), (max-width: 1199px) 30vw, 420px')
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
    fragment<-if(grepl('#',href,fixed=TRUE)) sub('^[^#]+','',href) else ''
    if(is.na(topic)) {
     routes <- page_migrations[[paste0('reference/',name,'.html')]]
     if(is.null(routes)) stop('Unresolved current-package autolink: ',href)
     key <- sub('^#','',fragment)
     position <- match(key, names(routes), nomatch=match('',names(routes)))
     destination <- routes[[position]]
    } else {
     destination <- paste0('reference/',topic,'.html',fragment)
    }
    xml2::xml_set_attr(a,'href',paste0(prefix,destination))
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
  # A deterministic local full-text index also supports Han substring search.
  xml2::write_html(doc,path)
 }
}
cat('Paired navigation and local search attached to all documentation pages.\n')
