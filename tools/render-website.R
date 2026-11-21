# Render articles and API reference bodies with the same pkgdown template.
# Article HTML is already converted from paired Markdown; no vignette is knitted.
Sys.setenv(CI = 'true') # init_site must use the local favicon, not its optional remote service.
options(sass.cache = FALSE)
page_index <- jsonlite::fromJSON('website/pages.json', simplifyVector = FALSE)
navigation <- jsonlite::fromJSON('website/navigation.json', simplifyVector = FALSE)
primary_groups <- Filter(function(group) !is.null(group$navbar), navigation)
arguments <- commandArgs(trailingOnly=TRUE)
option <- function(key, default) {
  value <- grep(paste0('^--',key,'='),arguments,value=TRUE)
  if(length(value)) sub('^[^=]+=','',value[[1]]) else default
}
part <- match.arg(option('part', 'articles'), c('init', 'articles'))
site <- option('site', '.site')
if (part == 'articles') {
  # Downlit resolves library() calls in the examples. Use this checkout's
  # metadata so highlighting works offline even before ggchord is installed.
  pkgload::load_all(".", export_all = FALSE, quiet = TRUE)
}
payload <- if (part == 'articles')
  jsonlite::fromJSON(option('payload', file.path(site, 'page-content.json')), simplifyVector = FALSE)
offset <- as.integer(option('offset','0'))
count <- as.integer(option('count','10000'))
languages <- arguments[!startsWith(arguments,'--')]
if(!length(languages)) languages <- c('en','zh')
stopifnot(all(languages %in% c('en','zh')))
for (lang in languages) {
  is_han <- lang == 'zh'
  Sys.setenv(LANGUAGE=if(is_han)'zh_CN' else 'en')
  keys <- vapply(primary_groups, function(group) group$navbar, '')
  components <- setNames(lapply(primary_groups, function(group) {
    page <- Filter(function(x) identical(x$slug, group$landing), page_index)[[1L]]
    list(text = group$title[[lang]], href = page$canonical_path)
  }), keys)
  cfg <- list(destination=file.path(site,lang),lang=if(is_han)'zh_CN' else 'en',
    navbar=list(structure = list(left = keys, right = c('github', 'lightswitch')),
      components=components))
  cfg$url <- paste0(yaml::read_yaml('_pkgdown.yml')$url, lang, '/')
  pkg <- pkgdown::as_pkgdown('.', override=cfg)
  pkg$meta$navbar$structure <- cfg$navbar$structure
  for(key in keys) pkg$meta$navbar$components[[key]] <- components[[key]]
  if (!file.exists(file.path(pkg$dst_path, 'lightswitch.js'))) pkgdown::init_site(pkg)
  if (part == 'init') {
    # build_home() also renders arbitrary root Markdown and the README. Those
    # pages have paired handbook sources; only the pkgdown utility pages apply.
    for (builder in c("build_citation_authors", "build_home_license", "build_404"))
      getFromNamespace(builder, "pkgdown")(pkg)
    next
  }
  selected <- payload[[lang]]
  if(length(selected)) selected <- selected[seq.int(offset+1L,min(length(selected),offset+count))]
  for(page in selected) {
    cat('Article:',lang,page$slug,'\n')
    path <- paste0(page$slug, '.html')
    dir.create(dirname(file.path(pkg$dst_path, path)), recursive = TRUE, showWarnings = FALSE)
    heading_id <- if (!is.null(page$heading_id))
      paste0(' id="', htmltools::htmlEscape(page$heading_id, attribute = TRUE), '"') else ''
    tryCatch(pkgdown::render_page(pkg,'home',list(pagetitle=page$title,
      index=paste0('<h1',heading_id,'>',htmltools::htmlEscape(page$title),'</h1>',page$body),
      sidebar=page$sidebar),path=path,quiet=TRUE),
      error=function(e) {cat('FAILED:',lang,page$slug,conditionMessage(e),'\n');stop(e)})
    cat('Completed article:',lang,page$slug,'\n')
    invisible(gc())
  }
}
cat('Rendered documentation with the shared pkgdown templates.\n')
