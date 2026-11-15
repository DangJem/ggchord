# Native pkgdown templates and reference with real example output.
# Article HTML is already converted from paired Markdown; no vignette is knitted.
Sys.setenv(CI = 'true') # init_site must use the local favicon, not its optional remote service.
options(sass.cache = FALSE)
page_index <- jsonlite::fromJSON('website/pages.json', simplifyVector = FALSE)
arguments <- commandArgs(trailingOnly=TRUE)
option <- function(key, default) {
  value <- grep(paste0('^--',key,'='),arguments,value=TRUE)
  if(length(value)) sub('^[^=]+=','',value[[1]]) else default
}
part <- match.arg(option('part','all'),c('all','reference','articles'))
# Article payloads are already HTML and never evaluate package code.
if (part != "articles") pkgload::load_all(quiet = TRUE)
payload <- if (part == 'reference') list(en = list(), zh = list()) else
  jsonlite::fromJSON('.site/page-content.json', simplifyVector = FALSE)
offset <- as.integer(option('offset','0'))
count <- as.integer(option('count','10000'))
languages <- arguments[!startsWith(arguments,'--')]
if(!length(languages)) languages <- c('en','zh')
stopifnot(all(languages %in% c('en','zh')))
reference_examples <- jsonlite::fromJSON('website/reference-examples.json')
if('zh' %in% languages && part!='articles') {
 status <- system2('python3',c('tools/reference-translations.py','--output','.site/reference-zh'))
 if(status != 0L) stop('Han reference translations are incomplete or stale')
}
for (lang in languages) {
  is_han <- lang == 'zh'
  Sys.setenv(LANGUAGE=if(is_han)'zh_CN' else 'en')
  nav_names <- if(is_han) c('首页','开始使用','绘图模块','通用技巧','数据与工作流','案例图库','API 参考') else
    c('Home','Get started','Plotting modules','Shared techniques','Data & workflows','Gallery','API reference')
  keys <- c('home','start','plots','techniques','data','gallery','reference')
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
  components$home$href <- 'index.html'
  components$start$menu <- page_menu(c('installation','modes','data'))
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
  if(is_han && part != "articles") {
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
  # Shared executable supplements cover topics with no runnable Rd example.
  # Keep code and native output together through pkgdown's ordinary evaluator.
  supplemental_topics <- names(reference_examples)
  if (part == "reference") {
    selected_files <- pkg$topics$file_in[seq.int(offset + 1L, min(nrow(pkg$topics), offset + count))]
    supplemental_topics <- intersect(supplemental_topics, tools::file_path_sans_ext(selected_files))
  }
  if (part != "articles") for (topic in supplemental_topics) {
    i <- match(paste0(topic, '.Rd'), pkg$topics$file_in)
    if (is.na(i)) stop('Unknown reference example topic: ', topic)
    input <- if (is_han) file.path('.site/reference-zh', pkg$topics$file_in[i]) else
      file.path('man', pkg$topics$file_in[i])
    text <- paste(readLines(input, warn = FALSE), collapse = "\n")
    at <- regexpr("\\examples{", text, fixed = TRUE)
    if (at[1L] > 0L) {
      chars <- strsplit(text, "", fixed = TRUE)[[1L]]
      cursor <- at[1L] + attr(at, "match.length")
      depth <- 1L
      while (depth > 0L) {
        if (chars[cursor] == "\\") cursor <- cursor + 1L else
          depth <- depth + (chars[cursor] == "{") - (chars[cursor] == "}")
        cursor <- cursor + 1L
      }
      text <- paste0(substr(text, 1L, at[1L] - 1L), substring(text, cursor))
    }
    code <- paste0("\\examples{\n", reference_examples[[topic]], "\n}\n")
    path <- tempfile(fileext = '.Rd')
    writeLines(c(text, code), path)
    pkg$topics$rd[[i]] <- getFromNamespace('rd_file', 'pkgdown')(path, pkg_path = '.')
    unlink(path)
  }
  if(part!='articles') {
  if (part != 'reference' || offset == 0L) {
  pkgdown::init_site(pkg)
  pkgdown::build_home(pkg, preview=FALSE)
  # build_home also converts development Markdown from the repository root.
  # Retain only package legal/author utilities; paired articles are rendered below.
  auxiliary <- c('index.html','authors.html','LICENSE.html','LICENSE-text.html','404.html')
  generated <- list.files(pkg$dst_path, '\\.html$', full.names=TRUE)
  unlink(generated[!basename(generated) %in% auxiliary])
  }
  reference_topics <- if (part == 'reference')
    pkg$topics$name[seq.int(offset + 1L, min(nrow(pkg$topics), offset + count))] else NULL
  # pkgdown evaluates under reference/, while source examples use the checkout's
  # relative input paths. Expose those inputs only during evaluation.
  input_link <- file.path(pkg$dst_path, 'reference', 'examples')
  dir.create(dirname(input_link), recursive = TRUE, showWarnings = FALSE)
  if (!file.symlink(normalizePath('examples'), input_link)) stop('Cannot expose local example inputs')
  tryCatch(pkgdown::build_reference(pkg,examples=TRUE,run_dont_run=FALSE,lazy=FALSE,preview=FALSE,
    topics = reference_topics),
    finally = unlink(input_link))
  }
  selected <- if(part=='reference') list() else payload[[lang]]
  if(length(selected)) selected <- selected[seq.int(offset+1L,min(length(selected),offset+count))]
  for(page in selected) {
    cat('Article:',lang,page$slug,'\n')
    tryCatch(pkgdown::render_page(pkg,'home',list(pagetitle=page$title,
      index=paste0('<h1>',htmltools::htmlEscape(page$title),'</h1>',page$body),
      sidebar=page$sidebar),path=paste0(page$slug,'.html'),quiet=TRUE),
      error=function(e) {cat('FAILED:',lang,page$slug,conditionMessage(e),'\n');stop(e)})
    cat('Completed article:',lang,page$slug,'\n')
    invisible(gc())
  }
}
cat('pkgdown templates and reference rendered for en and zh. Reference examples evaluated; no vignettes executed.\n')
