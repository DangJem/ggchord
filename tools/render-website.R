# Native pkgdown templates and complete reference, without executing examples.
# Article HTML is already converted from paired Markdown; no vignette is knitted.
Sys.setenv(CI = 'true') # init_site must use the local favicon, not its optional remote service.
options(sass.cache = FALSE)
pkgload::load_all(quiet = TRUE)
payload <- jsonlite::fromJSON('.site/page-content.json', simplifyVector = FALSE)
languages <- commandArgs(trailingOnly=TRUE)
if(!length(languages)) languages <- c('en','zh')
stopifnot(all(languages %in% c('en','zh')))
if('zh' %in% languages) {
 status <- system2('python3',c('tools/reference-translations.py','--output','.site/reference-zh'))
 if(status != 0L) stop('Han reference translations are incomplete or stale')
}
for (lang in languages) {
  is_han <- lang == 'zh'
  Sys.setenv(LANGUAGE=if(is_han)'zh_CN' else 'en')
  nav_names <- if(is_han) c('首页','开始使用','绘图章节','专题','函数参考') else
    c('Home','Get started','Plot chapters','Topics','Reference')
  keys <- c('home','start','plots','articles','reference')
  components <- setNames(lapply(seq_along(keys),function(i) list(text=nav_names[i])),keys)
  # The three task modules are shared by the top menu and the full site navigation.
  groups <- jsonlite::fromJSON('website/navigation.json', simplifyVector=FALSE)
  modules <- Filter(function(x) x$id %in% c('chord','collinear','plasmid'),groups)
  components$plots$menu <- lapply(modules,function(group) list(
    text=group$title[[lang]],href=paste0(group$pages[[1]],'.html')))
  topics <- Filter(function(x) x$id %in% c('data','techniques','cases'), groups)
  components$articles$menu <- c(list(list(
    text=if(is_han) '章节目录' else 'Chapter index', href='about.html')),
    lapply(topics, function(group) list(text=group$title[[lang]], href=paste0('about.html#path-',group$id))))
  cfg <- list(destination=paste0('.site/',lang),lang=if(is_han)'zh_CN' else 'en',
    navbar=list(components=components))
  # Localized category titles; topic bindings remain the exact source names.
  source <- yaml::read_yaml('_pkgdown.yml')
  cfg$url <- paste0(source$url,lang,'/')
  if(is_han) {
    titles <- c('创建图形','读取和整理数据','图层和标签','配色、图例和布局调整',
      '添加数据轨道','查找序列注释','投影和局部比较','保存与检查图形','方法、示例数据和通用工具')
    cfg$reference <- source$reference
    for(i in seq_along(cfg$reference)) cfg$reference[[i]]$title <- titles[i]
  }
  pkg <- pkgdown::as_pkgdown('.',override=cfg)
  # utils::modifyList does not replace the unnamed reference section list.
  # Assign the reviewed localized sections before the public reference builder.
  if(is_han) pkg$meta$reference <- cfg$reference
  pkg$meta$navbar$components$plots$menu <- components$plots$menu
  pkg$meta$navbar$components$articles <- components$articles
  if(is_han) {
    for(i in seq_len(nrow(pkg$topics))) {
      path <- file.path('.site/reference-zh',pkg$topics$file_in[i])
      raw <- tools::parse_Rd(path)
      # pkgdown topics require its classed Rd nodes, not tools::parse_Rd output.
      # Keep this one parser adapter explicit; generation uses public builders.
      reader <- getFromNamespace('rd_file','pkgdown')
      rd <- reader(path, pkg_path = '.')
      stopifnot(inherits(rd,'tag'), any(vapply(rd,inherits,logical(1),'tag_title')))
      pkg$topics$rd[[i]] <- rd
      title <- Filter(function(x) identical(attr(x,'Rd_tag'),'\\title'),raw)
      pkg$topics$title[i] <- paste(unlist(title),collapse='')
    }
  }
  pkgdown::init_site(pkg)
  pkgdown::build_home(pkg, preview=FALSE)
  # build_home also converts development Markdown from the repository root.
  # Retain only package legal/author utilities; paired articles are rendered below.
  auxiliary <- c('index.html','authors.html','LICENSE.html','LICENSE-text.html','404.html')
  generated <- list.files(pkg$dst_path, '\\.html$', full.names=TRUE)
  unlink(generated[!basename(generated) %in% auxiliary])
  pkgdown::build_reference(pkg,examples=FALSE,lazy=FALSE,preview=FALSE)
  for(page in payload[[lang]]) {
    cat('Article:',lang,page$slug,'\n')
    pkgdown::render_page(pkg,'home',list(pagetitle=page$title,
      index=paste0('<h1>',htmltools::htmlEscape(page$title),'</h1>',page$body),
      sidebar=page$sidebar),path=paste0(page$slug,'.html'),quiet=TRUE)
  }
}
cat('pkgdown templates and reference rendered for en and zh. No examples or vignettes executed.\n')
