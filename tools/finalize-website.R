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
navigation <- jsonlite::fromJSON('website/navigation.json', simplifyVector=FALSE)
page_index <- jsonlite::fromJSON('website/pages.json', simplifyVector=FALSE)
page_record <- function(key) {
  found <- Filter(function(x) identical(x$slug,key),page_index)
  if(length(found)) found[[1]] else NULL
}
figure_lessons <- jsonlite::fromJSON('website/figure-lessons.json')
module_motifs <- jsonlite::fromJSON('website/module-motifs.json')
course_stages <- jsonlite::fromJSON('website/course-stages.json', simplifyVector=FALSE)
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
    id <- if(title %in% c('绘图模块', 'Plotting modules')) 'site-visualizations'
      else if(title %in% c('通用技巧', 'Shared techniques')) 'site-techniques'
      else if(title %in% c('数据与工作流', 'Data & workflows')) 'site-data' else NULL
    if(is.null(id)) next
    old_id <- xml2::xml_attr(menu,'id')
    xml2::xml_set_attr(menu,'id',id)
    labelled <- xml2::xml_find_all(xml2::xml_parent(menu),'.//*[@aria-labelledby]')
    for(node in labelled) if(identical(xml2::xml_attr(node,'aria-labelledby'),old_id))
      xml2::xml_set_attr(node,'aria-labelledby',id)
  }
  main<-xml2::xml_find_first(doc,'//main')
  xml2::xml_set_attr(main,'tabindex','-1')
  slug <- sub('\\.html$','',basename(relative))
  if(slug=='LICENSE') {
    label <- if(lang=='zh') 'MIT 许可' else 'MIT License'
    xml2::xml_set_text(xml2::xml_find_first(main,'.//h1'),label)
    xml2::xml_set_text(xml2::xml_find_first(doc,'//title'),paste(label,'• ggchord'))
  }
  if(slug=='404' && lang=='zh') for(node in xml2::xml_find_all(main,'.//p'))
    xml2::xml_set_text(node,'未找到页面。请从首页、学习总览或搜索找回所需内容。')
  xml2::xml_set_attr(xml2::xml_find_first(doc,'//body'),'data-page',if(startsWith(relative,'reference/')) 'reference' else slug)
  if(slug=='index' && length(xml2::xml_find_all(main,'.//h1'))>1L) xml2::xml_remove(xml2::xml_find_first(main,'./h1'))
  layout <- xml2::xml_parent(main)
  xml2::xml_set_attr(layout,'class',paste(xml2::xml_attr(layout,'class'),'doc-layout'))
  prefix <- paste(rep('../',deep),collapse='')
  current_page <- if(startsWith(relative, 'reference/')) 'reference' else slug
  page_link <- function(key, label = NULL) {
    page <- Filter(function(x) identical(x$slug, key), page_index)[[1]]
    target <- page$canonical_path
    if(is.null(label)) label <- page$title[[lang]]
    paste0('<a ', if(key == current_page) 'aria-current="page" ' else '',
      'href="', prefix, target, '">', htmltools::htmlEscape(label), '</a>')
  }
  selected_groups <- Filter(function(group) current_page %in% group$pages, navigation)
  if(!length(selected_groups)) selected_groups <- Filter(function(group) group$id=='information',navigation)
  chapter_links <- paste(vapply(selected_groups, function(group) {
    stages <- course_stages[[group$id]]
    links <- ''
    for(key in unlist(group$pages)) {
      stage <- if(is.null(stages)) NULL else Filter(function(x) key %in% unlist(x$pages),stages)
      if(length(stage) && identical(key,stage[[1]]$pages[[1]])) links <- paste0(links,
        '<span class="learning-stage">',htmltools::htmlEscape(stage[[1]]$title[[lang]]),'</span>')
      links <- paste0(links,page_link(key))
    }
    paste0('<a class="course-home" href="',prefix,sub('[.]html$','',page_record(group$landing)$canonical_path),'.html">',
      if(lang=='zh')'课程首页 · ' else 'Course home · ',htmltools::htmlEscape(group$title[[lang]]),
      '</a><details class="chapter-section" open><summary>',htmltools::htmlEscape(group$title[[lang]]),
      '</summary>',links,'</details><a class="other-courses" href="',prefix,'learn.html">',
      if(lang=='zh')'浏览其他课程与专题 →' else 'Browse other courses and topics →','</a>')
  }, ''), collapse = '')
  nav_fragment <- function(html) xml2::xml_find_first(xml2::read_html(paste0('<html><body>',html,'</body></html>')),'//body/*')
  xml2::xml_add_sibling(main,nav_fragment(paste0('<nav class="chapter-nav" aria-label="',
    if(lang=='zh')'章节与参考' else 'Chapters and reference','"><p class="navigation-title">',
    if(lang=='zh')'章节与参考' else 'Chapters and reference','</p>',chapter_links,'</nav>')),.where='before')
  xml2::xml_add_child(main,nav_fragment(paste0('<details class="chapter-drawer"><summary>',
    if(lang=='zh')'章节与参考' else 'Chapters and reference','</summary><nav>',chapter_links,'</nav></details>')),.where=0L)

  current_group <- Filter(function(group) current_page %in% group$pages, navigation)
  group <- selected_groups[[1]]
  record <- page_record(current_page)
  module <- if(!is.null(record)) record$module else 'information'
  parent <- if(startsWith(relative,'reference/')) 'reference' else if(!is.null(record$parent)) record$parent else 'about'
  parent_record <- page_record(parent)
  xml2::xml_set_attr(xml2::xml_find_first(doc,'//body'),'data-module',module)
  crumb <- paste0('<nav class="breadcrumb" aria-label="',if(lang=='zh')'位置导航' else 'Breadcrumb','">',
    '<a href="',prefix,'index.html">',if(lang=='zh')'首页' else 'Home','</a><span aria-hidden="true">/</span>',
    '<a href="',prefix,parent_record$canonical_path,'">',htmltools::htmlEscape(parent_record$title[[lang]]),
    '</a><span aria-hidden="true">/</span><span aria-current="page">',
    htmltools::htmlEscape(xml2::xml_text(xml2::xml_find_first(main,'.//h1'))),'</span></nav>')
  if(current_page!='index') xml2::xml_add_sibling(xml2::xml_find_first(main,'.//h1'),nav_fragment(crumb),.where='before')
  if(length(current_group) && current_page != 'index' && !startsWith(relative, 'reference/')) {
    group <- current_group[[1]]
    group_title <- htmltools::htmlEscape(group$title[[lang]])
    sequence <- !is.null(record) && length(record$next_steps)>0L
    stage <- group_title
    if(sequence) {
      step <- match(current_page, group$pages)
      total <- length(group$pages)
      stage <- paste(group_title, '·', if(lang == 'zh') sprintf('第 %d 节，共 %d 节', step, total)
        else sprintf('Section %d of %d', step, total))
      prev <- if(step > 1) page_link(group$pages[[step - 1]], if(lang == 'zh') '← 上一节' else '← Previous section')
      next_link <- page_link(record$next_steps[[1]], if(lang == 'zh') '下一步 →' else 'Next step →')
      end <- page_link('learn', if(lang == 'zh') '学习总览' else 'Learning paths')
      xml2::xml_add_child(main, nav_fragment(paste0('<nav class="lesson-pagination" aria-label="',
        if(lang == 'zh') '章节进度' else 'Chapter navigation', '">', prev,
        next_link, '</nav>')))
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
  if(module %in% names(module_motifs)) {
    symbol <- nav_fragment(paste0('<span class="course-symbol"><svg viewBox="0 0 96 48" aria-hidden="true" focusable="false">',module_motifs[[module]],'</svg></span>'))
    xml2::xml_add_child(xml2::xml_find_first(main,'.//h1'),symbol,.where=0L)
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
  brand <- xml2::xml_find_first(doc,'//a[contains(@class,"navbar-brand")]')
  xml2::xml_add_child(brand,nav_fragment(paste0('<img src="',base,'logo.webp" width="36" height="36" alt="" aria-hidden="true">')),.where=0L)
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
  search<-fragment(paste0('<li class="nav-item local-search"><button type="button" class="search-trigger" id="site-search-trigger" aria-haspopup="dialog"><span aria-hidden="true">⌕</span><span class="search-trigger-label">',
   if(lang=='zh')'搜索任务、函数、参数' else 'Search tasks, functions, parameters','</span><kbd aria-hidden="true">⌘ K</kbd></button>',
   '<noscript><a href="',prefix,'task-index.html">',label,'</a></noscript></li>'))
  xml2::xml_add_child(target,search)
  body<-xml2::xml_find_first(doc,'//body')
  xml2::xml_set_attr(body,'data-kind',if(startsWith(relative,'reference/')) 'api' else if(!is.null(record)) record$kind else 'utility')
  metadata <- jsonlite::fromJSON('.site/build-metadata.json')
  xml2::xml_add_child(target,fragment(paste0('<li class="nav-item nav-text"><a class="nav-link" href="',prefix,'installation.html">',
    metadata$package_version,' · ',metadata$channel,'</a></li>')))
  xml2::xml_add_child(body,fragment(paste0('<dialog class="search-dialog" id="site-search-dialog" aria-labelledby="site-search-title"><div class="search-heading"><h2 id="site-search-title">',label,
   '</h2><button type="button" class="search-close" aria-label="',if(lang=='zh')'关闭搜索' else 'Close search','">×</button></div>',
   '<label class="search-field" for="site-search"><span class="visually-hidden">',label,'</span><input id="site-search" type="search" role="combobox" aria-autocomplete="list" aria-controls="site-results" aria-expanded="false" autocomplete="off" placeholder="',
   if(lang=='zh')'试试：酶切引线、link_branch、图例' else 'Try: primer labels, link_branch, guides','"></label><div class="search-filters" aria-label="',
   if(lang=='zh')'筛选结果类型' else 'Filter result types','"></div><p class="search-count" id="search-status" role="status"></p><ul id="site-results" class="search-results" role="listbox" aria-label="',
   if(lang=='zh')'搜索结果' else 'Search results','"></ul><div class="search-help">',
   if(lang=='zh')'↑ ↓ 选择 · Enter 打开 · Esc 关闭；也可按任务查找。' else '↑ ↓ select · Enter open · Esc close; or browse by task.',
   ' <a href="',prefix,'task-index.html">',if(lang=='zh')'任务索引' else 'Task index','</a></div></dialog>')))
  xml2::xml_set_attr(body,'data-lang',lang)
  xml2::xml_set_attr(body,'data-search',paste0(paste(rep('../',deep),collapse=''),'search.json'))
  if(startsWith(relative,'reference/')) {
   info<-if(lang=='zh') '本页完整说明已译为现代汉语；函数、参数、代码和原始数据标识保留原样。' else 'Complete canonical API text; task articles explain the workflows.'
   reviews <- jsonlite::fromJSON('website/teaching-reviews.json',simplifyVector=FALSE)
   matching <- Filter(function(name) identical(alias_topics[[name]],slug),names(reviews))
   lesson <- if(length(matching)) reviews[[matching[[1]]]] else NULL
   target_page <- if(!is.null(lesson)) paste0(lesson$page,'.html#',lesson$anchor) else 'learn.html'
   link<-paste0(paste(rep('../',deep),collapse=''),target_page)
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
  if(file.exists(file.path('.site/assets/examples',paste0(current_page,'-',lang,'.R'))) && !startsWith(relative,'reference/')) {
    file_url <- paste0(base,'examples/',current_page,'-',lang,'.R')
    action <- paste0('<div class="example-actions"><button type="button" class="copy-full-example" data-source="',file_url,'">',
      if(lang=='zh')'复制本页完整 R 代码' else 'Copy complete page R code','</button><a download href="',file_url,'">',
      if(lang=='zh')'下载 R 脚本' else 'Download R script','</a><small>',
      if(lang=='zh')'包含本页准备与全部步骤。文件示例须先下载输入；预览和保存调用保留。' else 'Includes preparation and every page step. File examples need the input bundle; preview and save calls are retained.','</small></div>')
    if(current_page=='index') xml2::xml_add_sibling(xml2::xml_find_first(main,'.//pre'),fragment(action),.where='before')
    else xml2::xml_add_sibling(xml2::xml_find_first(main,'.//h1'),fragment(action),.where='after')
  }
  for(pre in xml2::xml_find_all(main,'.//pre')) {
   parent<-xml2::xml_parent(pre)
   if(xml2::xml_name(parent)=='details')next
   is_r <- grepl('sourceCode|r-output|usage',paste(xml2::xml_attr(pre,'class'),xml2::xml_attr(parent,'class')))
   label <- if(is_r) if(lang=='zh')'R · 本步代码' else 'R · Step code' else if(lang=='zh')'输出 / 文本' else 'Output / text'
   details<-fragment(paste0('<details class="code-disclosure" data-code-kind="',if(is_r)'r' else 'output','" open><summary>',label,'</summary>',as.character(pre),'</details>'))
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
  for(output in xml2::xml_find_all(main,'.//*[contains(@class,"data-output")]')) {
    xml2::xml_set_attr(output,'tabindex','0')
    xml2::xml_set_attr(output,'role','region')
    xml2::xml_set_attr(output,'aria-label',if(lang=='zh')'R 数据输出，可横向滚动' else 'R data output, horizontally scrollable')
    for(cell in xml2::xml_find_all(output,'.//td')) if(grepl('^-?[0-9]+(\\.[0-9]+)?$',xml2::xml_text(cell))) xml2::xml_set_attr(cell,'class','numeric')
  }
  # Adjacent code/figure pairs form comparisons without changing code order.
  children <- xml2::xml_children(main)
  is_pair <- function(node) {
    next_node <- xml2::xml_find_first(node,'following-sibling::*[1]')
    identical(xml2::xml_name(node),'details') && grepl('code-disclosure',xml2::xml_attr(node,'class') %||% '') &&
      !inherits(next_node,'xml_missing') && startsWith(xml2::xml_attr(next_node,'data-figure') %||% '', 'compare')
  }
  skip <- 0L
  for(i in seq_along(children)) {
    if(i<=skip) next
    node <- children[[i]]
    if(!is_pair(node)) next
    pairs <- list(); cursor <- node
    while(!inherits(cursor,'xml_missing') && is_pair(cursor)) {
      fig <- xml2::xml_find_first(cursor,'following-sibling::*[1]')
      pairs[[length(pairs)+1L]] <- list(code=cursor,figure=fig)
      cursor <- xml2::xml_find_first(fig,'following-sibling::*[1]')
    }
    if(length(pairs)<2L) next
    skip <- i+2L*length(pairs)-1L
    group_node <- fragment('<div class="figure-comparison"></div>')
    xml2::xml_add_sibling(node,group_node,.where='before')
    group_node <- xml2::xml_find_first(node,'preceding-sibling::*[1]')
    for(pair in pairs) {
      xml2::xml_add_child(group_node,fragment('<div class="figure-panel"></div>'))
      panel <- xml2::xml_find_first(group_node,'./div[last()]')
      xml2::xml_add_child(panel,pair$code);xml2::xml_add_child(panel,pair$figure)
      xml2::xml_set_attr(xml2::xml_find_all(panel,'.//img'),'sizes','(max-width: 650px) calc(100vw - 48px), (max-width: 1199px) 42vw, 480px')
      xml2::xml_remove(pair$code);xml2::xml_remove(pair$figure)
    }
  }
  if(current_page=='parameter-atlas') {
    for(details in xml2::xml_find_all(main,'.//details')) {
      api <- sub('\\(\\)$','',xml2::xml_text(xml2::xml_find_first(details,'./summary/code')))
      for(row in xml2::xml_find_all(details,'.//table/tbody/tr')) {
        parameter <- xml2::xml_text(xml2::xml_find_first(row,'./td[1]'))
        xml2::xml_set_attr(row,'id',paste0('param-',api,'-',parameter))
      }
    }
  }
  for(figure in xml2::xml_find_all(main,'.//figure[@data-figure]')) {
    key <- xml2::xml_attr(figure,'data-figure')
    code_url <- if(current_page=='gallery') paste0(base,'documentation-examples.R') else
      paste0(prefix,if(grepl('#',figure_lessons[[key]],fixed=TRUE)) sub('#','.html#',figure_lessons[[key]],fixed=TRUE) else paste0(figure_lessons[[key]],'.html'))
    xml2::xml_set_attr(xml2::xml_find_first(figure,'./a'),'data-code-url',code_url)
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
