# Attach paired-language navigation and local progressive enhancements to pkgdown pages.
args<-commandArgs(trailingOnly=TRUE)
publication<-identical(args,'publication')
topic_files<-list.files('man','\\.Rd$',full.names=TRUE)
alias_topics<-unlist(lapply(topic_files,function(path) {
 rd<-tools::parse_Rd(path)
 aliases<-vapply(Filter(function(x) identical(attr(x,'Rd_tag'),'\\alias'),rd),
  function(x) paste(unlist(x),collapse=''),'')
 setNames(rep(sub('\\.Rd$','',basename(path)),length(aliases)),aliases)
}))
navigation <- jsonlite::fromJSON('website/navigation.json', simplifyVector=FALSE)
page_index <- jsonlite::fromJSON('website/pages.json', simplifyVector=FALSE)
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
  main<-xml2::xml_find_first(doc,'//main')
  xml2::xml_set_attr(main,'tabindex','-1')
  slug <- sub('\\.html$','',basename(relative))
  xml2::xml_set_attr(xml2::xml_find_first(doc,'//body'),'data-page',if(startsWith(relative,'reference/')) 'reference' else slug)
  layout <- xml2::xml_parent(main)
  xml2::xml_set_attr(layout,'class',paste(xml2::xml_attr(layout,'class'),'doc-layout'))
  prefix <- paste(rep('../',deep),collapse='')
  chapter_links <- paste(vapply(navigation,function(group) {
    links <- paste(vapply(group$pages,function(key) {
      page <- Filter(function(x) identical(x$slug,key),page_index)[[1]]
      target <- if(key=='reference') 'reference/index.html' else paste0(key,'.html')
      active <- key==slug || (key=='reference' && startsWith(relative,'reference/'))
      paste0('<a ',if(active)'aria-current="page" ' else '', 'href="',prefix,target,'">',
        htmltools::htmlEscape(page$title[[lang]]),'</a>')
    },''),collapse='')
    paste0('<section><p class="chapter-group">',group$title[[lang]],'</p>',links,'</section>')
  },''),collapse='')
  nav_fragment <- function(html) xml2::xml_find_first(xml2::read_html(paste0('<html><body>',html,'</body></html>')),'//body/*')
  xml2::xml_add_sibling(main,nav_fragment(paste0('<nav class="chapter-nav" aria-label="',
    if(lang=='zh')'文档章节' else 'Documentation chapters','">',chapter_links,'</nav>')),.where='before')
  xml2::xml_add_child(main,nav_fragment(paste0('<details class="chapter-drawer"><summary>',
    if(lang=='zh')'浏览文档' else 'Browse documentation','</summary><nav>',chapter_links,'</nav></details>')),.where=0L)
  # Native reference TOCs are populated later by JavaScript. Build the compact
  # outline from actual headings, rather than cloning their initially empty nav.
  outline_headings <- xml2::xml_find_all(main,'.//*[self::h2 or self::h3][@id]')
  if(length(outline_headings)) {
    outline_links <- paste(vapply(outline_headings,function(h) paste0('<li><a href="#',
      htmltools::htmlEscape(xml2::xml_attr(h,'id'),attribute=TRUE),'">',
      htmltools::htmlEscape(xml2::xml_text(h)),'</a></li>'),''),collapse='')
    toc_fragment <- nav_fragment(paste0('<details class="page-outline"><summary>',
      if(lang=='zh')'本页目录' else 'On this page','</summary><ul>',outline_links,'</ul></details>'))
    xml2::xml_add_sibling(xml2::xml_find_first(main,'.//h1'),toc_fragment,.where='after')
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
  language<-fragment(paste0('<a id="site-language" class="language-control" href="',partner,
   '" lang="',if(lang=='en')'zh-Hans' else 'en','">',if(lang=='en')'现代汉语' else 'English','</a>'))
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
  notice<-if(lang=='zh') '0.16.5 开发版文档 · 尚未正式发布' else '0.16.5 development documentation · unreleased'
  if(!publication) notice<-paste(notice,if(lang=='zh')'· 本地预览，尚未部署' else '· local preview, not deployed')
  xml2::xml_add_child(main,fragment(paste0('<div class="candidate-note">',notice,'</div>')),.where=0L)
  if(startsWith(relative,'reference/')) {
   info<-if(lang=='zh') '本页完整说明已译为现代汉语；函数、参数、代码和原始数据标识保留原样。' else 'Complete canonical API text; task articles explain the workflows.'
   link<-paste0(paste(rep('../',deep),collapse=''),'articles.html')
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
   if(grepl('/assets/(NEWS|ASSET-SOURCES)\\.html$',href))
    xml2::xml_set_attr(a,'href',sub('\\.html$','.md',href))
   if(grepl('^https?://rdrr\\.io/pkg/ggchord/man/[^#]+\\.html',href)) {
    name<-utils::URLdecode(sub('\\.html.*$','',sub('^.*/man/','',href)))
    topic<-if(name %in% names(alias_topics)) alias_topics[[name]] else
     if(name %in% sub('\\.Rd$','',basename(topic_files))) name else NA_character_
    if(is.na(topic)) stop('Unresolved current-package autolink: ',href)
    fragment<-if(grepl('#',href,fixed=TRUE)) sub('^[^#]+','',href) else ''
    xml2::xml_set_attr(a,'href',paste0(paste(rep('../',deep),collapse=''),'reference/',topic,'.html',fragment))
   }
  }
  if(lang=='zh') {
   translations<-c('值'='返回值','范例'='示例','该页面'='本页目录')
   for(h in xml2::xml_find_all(doc,'//h2')) {
    children<-xml2::xml_contents(h)
    for(t in children) if(xml2::xml_type(t)=='text') {
     text<-trimws(xml2::xml_text(t))
     if(text %in% names(translations)) xml2::xml_set_text(t,translations[[text]])
    }
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
