# Native pkgdown templates and complete reference, without executing examples.
# Article HTML is already converted from paired Markdown; no vignette is knitted.
Sys.setenv(CI = 'true') # init_site must use the local favicon, not its optional remote service.
options(sass.cache = FALSE)
pkgload::load_all(quiet = TRUE)
payload <- jsonlite::fromJSON('.site/page-content.json', simplifyVector = FALSE)
languages <- commandArgs(trailingOnly=TRUE)
if(!length(languages)) languages <- c('en','zh')
stopifnot(all(languages %in% c('en','zh')))
for (lang in languages) {
  is_han <- lang == 'zh'
  Sys.setenv(LANGUAGE=if(is_han)'zh_CN' else 'en')
  nav_names <- if(is_han) c('入门','函数参考','专题','图集','变更') else
    c('Get started','Reference','Articles','Gallery','Changelog')
  keys <- c('start','reference','articles','gallery','news')
  components <- setNames(lapply(seq_along(keys),function(i) list(text=nav_names[i])),keys)
  cfg <- list(destination=paste0('.site/',lang),lang=if(is_han)'zh_CN' else 'en',
    navbar=list(components=components))
  # Localized category titles; topic bindings remain the exact source names.
  source <- yaml::read_yaml('_pkgdown.yml')
  cfg$url <- paste0(source$url,lang,'/')
  if(is_han) {
    titles <- c('构建与契约','输入与校验','几何与标签','位置、标度、主题与图例',
      '定量轨道','序列识别','投影与比较工作流','输出与诊断','方法、教学数据与上游重导出')
    cfg$reference <- source$reference
    for(i in seq_along(cfg$reference)) cfg$reference[[i]]$title <- titles[i]
  }
  pkg <- pkgdown::as_pkgdown('.',override=cfg)
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
