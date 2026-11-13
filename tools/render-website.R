# Native pkgdown templates and complete reference, without executing examples.
# Article HTML is already converted from paired Markdown; no vignette is knitted.
Sys.setenv(CI = 'true') # init_site must use the local favicon, not its optional remote service.
options(sass.cache = FALSE)
pkgload::load_all(quiet = TRUE)
payload <- jsonlite::fromJSON('.site/page-content.json', simplifyVector = FALSE)
page_index <- jsonlite::fromJSON('website/pages.json', simplifyVector = FALSE)
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
  nav_names <- if(is_han) c('开始使用','绘图模块','通用技巧','数据与工作流','案例图库','API 参考') else
    c('Get started','Plotting modules','Shared techniques','Data & workflows','Gallery','API reference')
  keys <- c('start','plots','techniques','data','gallery','reference')
  components <- setNames(lapply(seq_along(keys),function(i) list(text=nav_names[i])),keys)
  # Every primary menu opens useful teaching directly; the project page remains
  # an auxiliary index rather than a required intermediate navigation step.
  groups <- jsonlite::fromJSON('website/navigation.json', simplifyVector=FALSE)
  page_menu <- function(slugs) lapply(slugs, function(slug) {
    page <- Filter(function(x) identical(x$slug, slug), page_index)[[1L]]
    list(text = page$title[[lang]], href = page$canonical_path)
  })
  modules <- Filter(function(x) x$id %in% c('multi-sequence','circular','dotplot','fusion'),groups)
  components$plots$menu <- lapply(modules,function(group) list(
    text=group$title[[lang]],href=paste0(group$landing,'.html')))
  shared <- Filter(function(x) x$id == 'shared-techniques', groups)[[1L]]
  data <- Filter(function(x) x$id == 'data', groups)[[1L]]
  components$techniques$menu <- page_menu(unlist(shared$pages))
  components$data$menu <- page_menu(unlist(data$pages))
  components$start$href <- 'learn.html'
  components$gallery$href <- 'gallery.html'
  components$reference$href <- 'reference/index.html'
  cfg <- list(destination=paste0('.site/',lang),lang=if(is_han)'zh_CN' else 'en',
    navbar=list(structure = list(left = keys, right = c('github', 'lightswitch')),
      components=components))
  # Localized category titles; topic bindings remain the exact source names.
  source <- yaml::read_yaml('_pkgdown.yml')
  cfg$url <- paste0(source$url,lang,'/')
  if(is_han) {
    titles <- c('创建图形','读取和整理数据','图层和标签','配色、图例和布局调整',
      '添加数据轨道','查找序列注释','投影和局部比较','静态面板组合','保存与检查图形','方法、示例数据和通用工具')
    stopifnot(length(titles) == length(source$reference))
    cfg$reference <- source$reference
    for(i in seq_along(cfg$reference)) cfg$reference[[i]]$title <- titles[i]
  }
  pkg <- pkgdown::as_pkgdown('.',override=cfg)
  # utils::modifyList does not replace the unnamed reference section list.
  # Assign the reviewed localized sections before the public reference builder.
  if(is_han) pkg$meta$reference <- cfg$reference
  pkg$meta$navbar$structure <- cfg$navbar$structure
  for(key in keys) pkg$meta$navbar$components[[key]] <- components[[key]]
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
