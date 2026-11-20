# Native pkgdown templates and reference with real example output.
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
part <- match.arg(option('part','all'),c('all','init','reference','articles'))
site <- option('site', '.site')
# Article payloads are already HTML and never evaluate package code.
if (part %in% c("all", "reference")) pkgload::load_all(quiet = TRUE)
payload <- if (part %in% c('reference','init')) list(en = list(), zh = list()) else
  jsonlite::fromJSON(option('payload', file.path(site, 'page-content.json')), simplifyVector = FALSE)
offset <- as.integer(option('offset','0'))
count <- as.integer(option('count','10000'))
languages <- arguments[!startsWith(arguments,'--')]
if(!length(languages)) languages <- c('en','zh')
stopifnot(all(languages %in% c('en','zh')))
reference_examples <- if (part %in% c('all', 'reference')) jsonlite::fromJSON('website/reference-examples.json')
translated_reference <- file.path(site, 'reference-zh')
if('zh' %in% languages && part!='articles') {
 status <- system2('python3',c('tools/reference-translations.py','--output',shQuote(translated_reference)))
 if(status != 0L) stop('Han reference translations are incomplete or stale')
}
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
  # Localized category titles; topic bindings remain the exact source names.
  source <- yaml::read_yaml('_pkgdown.yml')
  if (part != 'articles') for (section in source$reference) {
    original <- paste(section$title, section$desc, sep = "\n")
    if (!identical(digest::digest(original, algo = "sha256", serialize = FALSE),
                   section$source_sha256) || !isTRUE(nzchar(section$title_zh)) ||
        !isTRUE(nzchar(section$desc_zh))) {
      stop("Missing or stale reference section translation: ", section$id)
    }
  }
  cfg$url <- paste0(source$url,lang,'/')
  if(is_han) {
    cfg$reference <- lapply(source$reference, function(section) {
      section$title <- section$title_zh
      section$desc <- section$desc_zh
      section
    })
  }
  pkg <- pkgdown::as_pkgdown('.',override=cfg)
  # utils::modifyList does not replace the unnamed reference section list.
  # Assign the reviewed localized sections before the public reference builder.
  if(is_han) pkg$meta$reference <- cfg$reference
  pkg$meta$navbar$structure <- cfg$navbar$structure
  for(key in keys) pkg$meta$navbar$components[[key]] <- components[[key]]
  if (part == 'articles' && !file.exists(file.path(pkg$dst_path, 'lightswitch.js'))) {
    pkgdown::init_site(pkg)
  }
  if(is_han && part != "articles") {
    for(i in seq_len(nrow(pkg$topics))) {
      path <- file.path(translated_reference,pkg$topics$file_in[i])
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
  if (part %in% c("all","reference")) for (topic in supplemental_topics) {
    i <- match(paste0(topic, '.Rd'), pkg$topics$file_in)
    if (is.na(i)) stop('Unknown reference example topic: ', topic)
    input <- if (is_han) file.path(translated_reference, pkg$topics$file_in[i]) else
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
  # Reference examples are copied independently of the console used by pkgdown.
  if (part %in% c("all", "reference")) {
    indices <- if (part == "reference") match(selected_files, pkg$topics$file_in) else seq_len(nrow(pkg$topics))
    for (i in indices) {
      rd <- pkg$topics$rd[[i]]
      j <- which(vapply(rd, inherits, logical(1), "tag_examples"))
      if (length(j) != 1L) stop("Missing reference example: ", pkg$topics$name[i])
      example <- rd[[j]]
      if (!any(trimws(unlist(example)) == "library(ggchord)")) {
        first <- structure("\nlibrary(ggchord)\n", class = "RCODE")
        rd[[j]] <- structure(c(list(first), unclass(example)), class = class(example))
        pkg$topics$rd[[i]] <- rd
      }
    }
  }
  if(part!='articles') {
  if (part %in% c('all','init')) {
  pkgdown::init_site(pkg)
  pkgdown::build_home(pkg, preview=FALSE)
  # build_home also converts development Markdown from the repository root.
  # Retain only package legal/author utilities; paired articles are rendered below.
  auxiliary <- c('index.html','authors.html','LICENSE.html','LICENSE-text.html','404.html')
  generated <- list.files(pkg$dst_path, '\\.html$', full.names=TRUE)
  unlink(generated[!basename(generated) %in% auxiliary])
  }
  if (part != 'init') {
  reference_topics <- if (part == 'reference')
    pkg$topics$name[seq.int(offset + 1L, min(nrow(pkg$topics), offset + count))] else NULL
  pkgdown::build_reference(pkg, examples = TRUE, run_dont_run = FALSE,
                           lazy = FALSE, preview = FALSE, topics = reference_topics)
  } else pkgdown::build_reference(pkg, examples=FALSE, topics=character(), preview=FALSE)
  }
  selected <- if(part %in% c('reference','init')) list() else payload[[lang]]
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
